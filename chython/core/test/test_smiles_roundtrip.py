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

`test_smiles_write_differential.py` puts the writer against oracles OUTSIDE the tree, which is the
stronger test of a convention -- a writer and a reader that share a misunderstanding agree with each
other and with nobody else -- so this file does not replace it.  It adds the one thing an external
oracle cannot give: a CLOSED loop, where the only two things in the circuit are ours, so a
disagreement is provably a bug in one of them and not a convention difference with somebody's
third-party parser.

WHAT THE LOOP CAN AND CANNOT PROVE.  `write_smiles(read_smiles(write_smiles(m, spec)), spec) ==
write_smiles(m, spec)` is a FIXED POINT claim about strings.  It catches a writer emitting something
the reader misreads and a reader building something the writer disagrees with.  It cannot catch a
misunderstanding both share, and it is not chemistry: two different strings can be the same
molecule, so a failure is only interesting once something INDEPENDENT says the molecule moved.  That
arbiter here is InChI, from `_inchi.pxi` -- a different canonicaliser with its own stereo layer,
which is exactly what is needed to tell "my canonical order picked the mirror labelling" from "the
round trip inverted a centre".  The stereo tests below use it and the constitution tests do not,
because the constitution corpora reach a fixed point outright and there is nothing to arbitrate.

AND THE INVARIANT CAN LIE IN ONE SHAPE, which is why `unknown_h` vetoes rather than merely explains:
the `A` dialect on a molecule whose hydrogen count nobody stated is a string fixed point AND a
changed molecule at the same time.  The last test in this file is that counterexample, worked in
full; the rule it establishes is that a fixed point is evidence and `unknown_h` is the veto.

MEASURED, on 2026-09-03, and every number below is asserted rather than described:

* `first_5K.smi`, 4999 public NCI records, four specs (`''`, `A`, `i`, `!s`): **4999/4999** a fixed
  point, zero outputs the reader refuses.  This corpus carries no stereo, so it is a pure
  constitution loop -- and it is the one that would catch a ring-closure digit or a bracket rule
  that only the writer believes in.
* `test/arenes.sdf` + `test/heterocycles_charges.smi`, 95 records: **41/41** of the records with a
  fully STATED hydrogen count are a fixed point.  The other 54 carry an atom whose count nobody
  stated, which SMILES has no way to spell, so the reader supplies a real count, the features move,
  the canonical order moves and the string moves with it.  That story is checked and not assumed:
  in STORED order, where the traversal does not consult the features at all, **94/95** are a fixed
  point.  Same molecules, same writer, and the difference is the ORDER and nothing else.
* `test/stereo.sdf`, 294 bridged records carrying 974 signs: **293/294** a fixed point as a string,
  and **293/293** of the stated-hydrogen records identical as a MOLECULE.  The one string that moves
  is this corpus's single unstated-hydrogen record; no record comes back with a byte-identical InChI
  and a different string, which is what a mirror automorphism surviving the search would look like.
* The NIBR PubChem examples, 340 bridged records emitting 1611 `/` and `\\` tokens: **340/340** a
  fixed point.  Double-bond geometry survives the loop outright, with no arbitration needed.

The corpora, the bridges and the rule about where a corpus may be read from all come from the
differential file; its fixtures are imported here rather than rebuilt, so there is one
definition of each corpus and one place to change when one of them moves.
"""
from pytest import mark

from chython.core import MoleculeContainer
from chython.core._core import molecule_to_inchi, read_smiles, smw_traversal, write_smiles
# imported for the side effect that pytest registers them as fixtures in THIS module too.  Also
# `bridge` and `bridge_stereo`, which are plain functions and are called directly.
#
# `oracle` is in the list because the four corpus fixtures REQUEST it: the corpora are read by an
# out-of-tree chython 2 co-process now, and a fixture that is not visible in the module that uses the
# ones depending on it is an ERROR at setup rather than a skip.  Leaving it out takes these ten tests
# out of the suite in a way that reads like an infrastructure problem instead of a missing import.
from chython.core.test.test_smiles_write_differential import (arenes, bridge, bridge_stereo,  # noqa
                                                              cis_trans, nci, oracle, stereo_sdf)


# `!s` is in the sweep although the corpus has no stereo: the spec still selects a different seed for
# the canonical order, so it is a different traversal reaching the same answer, and that is worth a
# column.  `r` is not here and never will be -- it raises.
SPECS = ['', 'A', 'i', '!s']


def loop(mol, spec=''):
    """`(out, out_again)`. The reader is given OUR string, so a refusal is our bug, not bad input."""
    out = write_smiles(mol, spec)
    return out, write_smiles(read_smiles(out), spec)


# ------------------------------------------------------------------------------------------------
# The constitution loop.  No stereo anywhere in this half, by measurement -- see the differential
# file's header -- so nothing here can pass or fail on a configuration.
@mark.parametrize('spec', SPECS)
def test_every_NCI_string_reads_back_to_ITSELF(nci, spec):
    """4999 of 4999, on all four specs, with the CLOSED loop: our reader in, our writer out.

    The source molecule is `read_smiles(text)` and not the V2 bridge, and that is the whole point of
    the word closed -- see the next test for what happens when one end of the circuit is somebody
    else's, and why the difference is not a defect in either.
    """
    moved = []
    refused = []
    total = 0
    for text, source in nci:
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


def test_with_CHYTHON_2_at_one_end_the_only_movers_are_the_UNSTATED_H_records(nci):
    """The open loop, and the exact price of SMILES having no spelling for "nobody said".

    Feed the writer a molecule that came from V2's parser instead of ours and two of 4999 records
    stop being a fixed point.  Both are reported by the writer in `unknown_h`, both are hypervalent
    nonsense (a ferrocene with nine bonds to iron, `C1=O=O1`), and the InChI of the one InChI will
    accept says `C11H46O3` before and `C11H16O3` after -- H_UNKNOWN is 15, twice over, so the phantom
    thirty hydrogens are the sentinel being read as a count by a tool that has no sentinel.  The V3
    reader never produces such an atom, which is exactly why the closed loop above does not see this at
    all.  Asserting `== 2` rather than "few" because a third would mean something new is unstated, and
    zero would mean the bridge stopped carrying the unstated state.

    A THIRD RECORD WAS A MOVER UNTIL THE CANONICAL ORDER BECAME PER COMPONENT: an aluminium urea
    complex with sulfate and `I[I]I`, where the unstated count moved the whole record's labelling and
    now moves only its own component's.  Confining the contamination to the component that has it is
    the property, not the count.
    """
    movers = []
    for text, source in nci:
        mol = bridge(source)
        out, again = loop(mol)
        if out != again:
            movers.append((text, out, again, smw_traversal(mol)['unknown_h']))
    assert len(movers) == 2, 'measured 2 movers on 2026-09-10, now %d: %r' % (
        len(movers), [m[0] for m in movers])
    explained = [m for m in movers if m[3]]
    assert len(explained) == 2, 'a mover the writer does NOT report as unstated: %r' % (
        [m[:3] for m in movers if not m[3]],)

    # and the same molecules in STORED order do not move at all, which localises the cause to the
    # canonical ORDER reacting to a hydrogen count rather than to anything about the tokens
    stored = [text for text, source in nci
              if loop(bridge(source), 'i')[0] != loop(bridge(source), 'i')[1]]
    assert not stored, 'stored order moved on %d records: %r' % (len(stored), stored[:3])


def test_the_aromatic_dialect_survives_its_own_UPPERCASE_brackets(arenes):
    """`A` writes `[CH]:1:[CH]` and not `c1cc`, deliberately, and the reader must accept that.

    `C:C` is outside OpenSMILES and no rule says what hydrogen count a bare aromatic upper-case atom
    implies, so the `A` dialect brackets the atom and states the count -- see
    `test_smiles_write_aromatic.py`.  The whole argument for that choice is that it is unambiguous,
    which is a claim about a READER, and until there was one the claim was untested.
    """
    refused = []
    bracketed = 0
    for text, source in arenes:
        out = write_smiles(bridge(source), 'A')
        if ':' in out:
            bracketed += 1
        try:
            read_smiles(out)
        except Exception as e:
            refused.append((out, '%s: %s' % (type(e).__name__, e)))
    # not every record in this corpus is STORED aromatic -- V2 hands some of them over kekulised, and
    # the `A` dialect writes what the arena holds rather than re-perceiving -- so the claim is about
    # the corpus and not about each record
    assert bracketed > 50, 'only %d of %d records produced an aromatic bond token, so this test was ' \
                           'not exercising the dialect' % (bracketed, len(arenes))
    assert not refused, '%d aromatic-dialect strings our reader refused: %r' % (len(refused),
                                                                                refused[:3])


def test_where_the_hydrogen_count_is_UNSTATED_the_ORDER_moves_and_not_the_molecule(arenes):
    """The 54 arene records that are not a fixed point, explained rather than excused.

    SMILES cannot spell "nobody stated a count", so the reader has to supply one; the count feeds the
    features, the features feed the canonical order, and the order picks a different start atom.  The
    check is that STORED order -- which asks the features nothing -- is a fixed point for all but the
    one record whose H count changes a ring's kekulisation, and that every canonical-order casualty
    is a record the writer ALREADY reported as unstated.  An unexplained casualty fails here.
    """
    unexplained = []
    stated = moved_stated = 0
    for text, source in arenes:
        mol = bridge(source)
        unknown = smw_traversal(mol)['unknown_h']
        out, again = loop(mol)
        if not unknown:
            stated += 1
            if out != again:
                moved_stated += 1
                unexplained.append((text, out, again))
    assert stated > 30, 'too few fully-stated records to mean anything: %d' % stated
    assert not moved_stated, ('%d of %d records with a stated hydrogen count moved: %r'
                              % (moved_stated, stated, unexplained[:3]))

    moved_stored = [t for t, s in arenes if loop(bridge(s), 'i')[0] != loop(bridge(s), 'i')[1]]
    assert len(moved_stored) <= 1, ('stored order should barely notice the unstated counts, and %d '
                                    'records moved: %r' % (len(moved_stored), moved_stored[:3]))


# ------------------------------------------------------------------------------------------------
# The configuration loop, and the half where the string moving is not the same question as the
# molecule moving.
def stereo_mols(records):
    """Bridged molecules from V2 records, dropping the ones the bridge itself cannot carry."""
    out = []
    for source in records:
        source = source[1] if isinstance(source, tuple) else source
        mol, carried, unset, failures = bridge_stereo(source)
        if not failures:
            out.append(mol)
    return out


def test_the_stereo_corpus_comes_back_as_the_SAME_MOLECULE(stereo_sdf):
    """293 of 293 stated-hydrogen records keep their InChI, stereo layer included.

    This is the assertion the whole file exists for.  It does not care which of an enantiomer's two
    labellings the canonical order chose, only that the compound that went in is the compound that
    came out -- and InChI, being nobody's SMILES, is entitled to that opinion.
    """
    mols = stereo_mols(stereo_sdf)
    assert len(mols) > 280, 'the bridge lost too many records: %d' % len(mols)
    assert sum(write_smiles(m).count('@') for m in mols) > 900, 'the corpus lost its signs'

    changed = []
    unstated = 0
    for mol in mols:
        out, again = loop(mol)
        if out == again:
            continue
        if smw_traversal(mol)['unknown_h']:
            unstated += 1                 # SMILES cannot carry it; the differential file measures it
            continue
        if molecule_to_inchi(mol) != molecule_to_inchi(read_smiles(out)):
            changed.append((out, again))
    assert unstated == 1, 'this corpus had exactly one unstated-hydrogen record, now %d' % unstated
    assert not changed, ('%d records came back as a DIFFERENT compound: %r'
                         % (len(changed), changed[:3]))


def test_every_string_is_a_FIXED_POINT_except_the_one_SMILES_cannot_state(stereo_sdf):
    """293 fixed points, 0 mirror-labelled, one unstated-hydrogen record -- pinned.

    The mirror column is EMPTY because the canonical search's leaf certificate carries a parity tail
    and its orbit prune refines by parity; 293 records are string fixed points.  Both numbers are
    pinned in both directions -- a rise means a mirror-labelled record came back, and a fall means a
    record stopped surviving its own round trip -- which makes this the second independent detector of
    the mirror-automorphism behaviour.

    Independent is the word that earns it a place next to the differential file.  That one sweeps
    creation orders explicitly; this one never permutes anything, because reading our own output IS a
    re-presentation in a different order -- the reader assigns slots by string position, which has no
    reason to agree with the order the writer was handed.  A fix that only satisfied a shuffle would
    still fail here.

    The remaining 1 is not stereo at all: it is the record whose hydrogen count nobody stated, which
    SMILES has no way to carry and which the last test in this file dissects on a benzene.
    """
    mols = stereo_mols(stereo_sdf)
    mirror = same = unstated = 0
    for mol in mols:
        out, again = loop(mol)
        if out == again:
            same += 1
        elif smw_traversal(mol)['unknown_h']:
            unstated += 1
        elif molecule_to_inchi(mol) == molecule_to_inchi(read_smiles(out)):
            mirror += 1
    assert (same, mirror, unstated) == (293, 0, 1), (
        'measured 293 fixed / 0 mirror-labelled / 1 unstated on 2026-09-03, now %r'
        % ((same, mirror, unstated),))


def test_every_cis_trans_record_is_a_fixed_point_with_nothing_to_arbitrate(cis_trans):
    """340 of 340, and 1611 direction tokens, so the `/` and `\\` half needs no InChI at all.

    Worth stating as its own test rather than folding into the tetrahedral one: a mirror automorphism
    is a TETRAHEDRAL story -- a reflection maps a stereocentre to its opposite and a double bond to
    itself -- so double-bond geometry has no reason to oscillate, and this test is the evidence for
    that reasoning rather than the reasoning on its own.
    """
    mols = stereo_mols(cis_trans)
    assert len(mols) > 300, 'the bridge lost too many records: %d' % len(mols)
    tokens = sum(write_smiles(m).count('/') + write_smiles(m).count('\\') for m in mols)
    assert tokens > 1500, 'the corpus lost its direction tokens: %d' % tokens

    moved = [(out, again) for out, again in (loop(m) for m in mols) if out != again]
    assert not moved, '%d of %d cis/trans records moved: %r' % (len(moved), len(mols), moved[:3])


# ------------------------------------------------------------------------------------------------
# WHERE THE FIXED POINT IS NOT ENOUGH.  The one shape in which this whole file's invariant lies.
def test_a_string_fixed_point_does_NOT_prove_the_molecule_survived():
    """The `A` dialect closes the loop on a benzene whose hydrogen count nobody stated -- and the
    molecule changes anyway.  This is the counterexample that justifies every `unknown_h` exclusion
    above.

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
