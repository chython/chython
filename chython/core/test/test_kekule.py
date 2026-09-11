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
"""`kekule()`: an aromatic edge set becomes bond orders 1 and 2.

Every molecule here is built through the container's own builder rather than parsed, because no
reader calls this -- a reader stores what its input said, and kekulisation is a deliberate operation
on top -- and because a builder fixture states exactly the facts a case is about.  That reason
outlived the one it was written with, which was that the SMILES reader did not exist yet.  The one
place a string is the better fixture is `test_the_six_vendor_spellings_end_to_end`, and it says why.

THE ARENA NOW STORES ORDER 4, and the tripwire that anticipated it has been spent: it was
`test_the_order_four_gate_is_still_closed`, it failed on the merge exactly as its comment said it
would, and it is replaced below by `test_the_no_argument_form_kekulises_the_stored_aromatic_bonds`.
The explicit edge set stays the dominant spelling here because most of these cases are about the
classifier's verdict on a stated subgraph, and because passing a SUBSET remains a real caller need
(a vendor file's per-bond aromatic mark).  `aromatic_bonds=None` now means the stored set rather
than nothing, which is `arom_stored_bonds`.
"""
from pytest import raises

from chython.core import MoleculeContainer
from chython.core._core import (AromaticKekulizeError, KekuleResult, kekule, kekule_classify,
                                kekule_copy)


def repairs(result):
    """The log without the `kekule:kekulized` line every system that comes out writes.

    That line is the WORK -- `INFO`, one per system, and the reason `mol.log` is not empty after a
    kekulisation that repaired nothing.  Every case below asks what its input NEEDED repairing, which
    is a different question, so it is filtered here rather than spelled into 39 assertions.  The record
    itself is pinned by `test_every_system_that_comes_out_says_so` and by
    `core/test/test_container_log.py`.
    """
    return [x for x in result.log if x.rule != 'kekule:kekulized']


def spread(result):
    """`(changed, repairs, unresolved)` from a `KekuleResult`, for the tests that read all three."""
    return result.changed, repairs(result), result.unresolved


def clean(result):
    """Something moved, and there was nothing to repair."""
    return result.changed and not repairs(result) and result.unresolved == []


def quiet(result):
    """Nothing moved at all -- so not even the kekulisation line, which only a rewritten system gets."""
    return not result.changed and result.log == [] and result.unresolved == []


def build(elements, aromatic, extra=(), charges=None, radicals=(), stated_h=None):
    """A molecule plus the aromatic edge set, in stable ids.

    `elements` is a sequence of symbols indexed 0..n-1 for the rest of the spec to refer to;
    `aromatic` and `extra` are `(i, j)` and `(i, j, order)` on those indices.  Returns
    `(mol, ids, aromatic_bonds, stated_h)` with everything already translated into stable ids,
    since that is what `kekule` takes.
    """
    mol = MoleculeContainer()
    charges = charges or {}
    ids = []
    for i, element in enumerate(elements):
        ids.append(mol.add_atom(element, charge=charges.get(i, 0), radical=i in radicals))
    bonds = []
    for i, j in aromatic:
        mol.add_bond(ids[i], ids[j], 1)
        bonds.append((ids[i], ids[j]))
    for i, j, order in extra:
        mol.add_bond(ids[i], ids[j], order)
    h = None if stated_h is None else {ids[i]: v for i, v in stated_h.items()}
    return mol, ids, bonds, h


def cycle(n, offset=0):
    return [(offset + i, offset + (i + 1) % n) for i in range(n)]


def orders(mol, bonds):
    return [mol.order_of(a, b) for a, b in bonds]


def alternating(mol, bonds):
    """Every atom of the system carries exactly one double bond.

    Checked as a property rather than against a specific bond pattern: which Kekule form comes
    out is a free choice between equivalent answers, and pinning one of them would make this test
    fail on a search reordering that is not a defect.
    """
    counts = {}
    for a, b in bonds:
        order = mol.order_of(a, b)
        assert order in (1, 2), (a, b, order)
        if order == 2:
            counts[a] = counts.get(a, 0) + 1
            counts[b] = counts.get(b, 0) + 1
    return counts


# --- the systems that must come out fully alternating

def test_benzene():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert log == []
    assert unresolved == []
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_naphthalene():
    aromatic = cycle(6) + [(4, 6), (6, 7), (7, 8), (8, 9), (9, 3)]
    mol, ids, bonds, _ = build('C' * 10, aromatic)
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_azulene():
    # a 5-ring and a 7-ring sharing an edge: the case a ring-by-ring kekuliser gets wrong,
    # because neither ring alone has an alternating assignment consistent with the other
    aromatic = cycle(5) + [(0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4)]
    mol, ids, bonds, _ = build('C' * 10, aromatic)
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_pyridine_unstated_hydrogen():
    # five must-match carbons is an odd count, so the may-match N is forced in and the answer is
    # pyridine without anything having to prefer it
    mol, ids, bonds, _ = build(['C'] * 5 + ['N'], cycle(6))
    # 'may' going IN, and the assertion belongs here rather than after the call: the free choice is
    # what the search is being given, and it is spent by being used
    assert kekule_classify(mol, bonds)[ids[5]] == 'may'
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}
    # AND 'must' COMING OUT, because the count is no longer unstated: the ring gave the nitrogen a
    # double bond, `kekule`'s hydrogen heal wrote the 0 that follows, and an atom with a stated 0 is
    # not a free choice any more.  A 'may' here would depend on the count staying open, which is
    # exactly what the heal closes, and the narrower contract is what stops a second call landing on
    # the other reading.
    assert mol.implicit_h_of(ids[5]) == 0
    assert kekule_classify(mol, bonds)[ids[5]] == 'must'


def test_pyridinium():
    mol, ids, bonds, _ = build(['C'] * 5 + ['N'], cycle(6), charges={5: 1})
    assert kekule_classify(mol, bonds)[ids[5]] == 'may'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_n_methylpyridinium_n_is_must():
    # charge +1 with three neighbours is pyridinium proper: it takes a ring double bond
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'C'], cycle(6), extra=[(5, 6, 1)],
                               charges={5: 1})
    assert kekule_classify(mol, bonds)[ids[5]] == 'must'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_pyrylium():
    mol, ids, bonds, _ = build(['C'] * 5 + ['O'], cycle(6), charges={5: 1})
    assert kekule_classify(mol, bonds)[ids[5]] == 'must'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids}


# --- the systems with a lone-pair donor, where one atom must stay single-bonded

def test_pyrrole_stated_nh():
    mol, ids, bonds, h = build(['C'] * 4 + ['N'], cycle(5), stated_h={4: 1})
    assert kekule_classify(mol, bonds, h)[ids[4]] == 'must_not'
    changed, log, unresolved = spread(kekule(mol, bonds, h))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_pyrrole_unstated_hydrogen():
    # four must-match carbons pair off among themselves, so the may-match N ends up unmatched and
    # the answer is pyrrole -- again without a preference rule
    mol, ids, bonds, _ = build(['C'] * 4 + ['N'], cycle(5))
    assert kekule_classify(mol, bonds)[ids[4]] == 'may'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_n_substituted_pyrrole():
    mol, ids, bonds, _ = build(['C'] * 4 + ['N', 'C'], cycle(5), extra=[(4, 5, 1)])
    assert kekule_classify(mol, bonds)[ids[4]] == 'must_not'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_furan():
    mol, ids, bonds, _ = build(['C'] * 4 + ['O'], cycle(5))
    assert kekule_classify(mol, bonds)[ids[4]] == 'must_not'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_thiophene():
    mol, ids, bonds, _ = build(['C'] * 4 + ['S'], cycle(5))
    assert kekule_classify(mol, bonds)[ids[4]] == 'must_not'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_pyrrolide():
    mol, ids, bonds, _ = build(['C'] * 4 + ['N'], cycle(5), charges={4: -1})
    assert kekule_classify(mol, bonds)[ids[4]] == 'must_not'
    assert clean(kekule(mol, bonds))


def test_cyclopentadienyl_anion():
    mol, ids, bonds, _ = build('C' * 5, cycle(5), charges={0: -1})
    assert kekule_classify(mol, bonds)[ids[0]] == 'may'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[1:]}


# --- an exocyclic double bond saturates its ring atom (MDL fixture 3)

def test_para_quinone_written_aromatic():
    mol, ids, bonds, _ = build(['C'] * 6 + ['O', 'O'], cycle(6),
                               extra=[(0, 6, 2), (3, 7, 2)])
    classes = kekule_classify(mol, bonds)
    assert classes[ids[0]] == 'must_not'
    assert classes[ids[3]] == 'must_not'
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[1:3] + ids[4:6]}


def test_pyridone_written_aromatic():
    # 2-pyridone spelled with an aromatic ring and an exocyclic carbonyl: the N donates its lone
    # pair and the carbonyl carbon is already saturated, so four carbons remain -- an even count
    aromatic = cycle(6)
    mol, ids, bonds, h = build(['N'] + ['C'] * 5 + ['O'], aromatic, extra=[(1, 6, 2)],
                               stated_h={0: 1})
    classes = kekule_classify(mol, bonds, h)
    assert classes[ids[0]] == 'must_not'
    assert classes[ids[1]] == 'must_not'
    assert clean(kekule(mol, bonds, h))
    assert alternating(mol, bonds) == {i: 1 for i in ids[2:6]}


# --- N-oxides, and the two gates that decide how hard this file tries to read one
#
# Every SMILES named in this section is a spelling files carry, and all six are read without one
# pattern per shape: the three rules below do the whole job.
#
# The rules being exercised, in the order the code applies them:
#
# * `n(=O)` and `n(=N)` on a NEUTRAL ring nitrogen are always rewritten charge-separated.  A neutral
#   aromatic N with two ring bonds has spent its lone pair on the ring, so the pi bond is not a
#   spelling anyone can honour.  The rewrite conserves total charge.
# * a BALANCED pair of mis-spellings in one aromatic system -- one nitrogen written cationic, one
#   written neutral-with-an-anion -- is taken together, because taken together it conserves charge.
# * a single mis-spelling in either direction MOVES the total charge, so it is offered only to a
#   system that has no Kekule form as written.  That gate, and nothing about ring size, is what
#   keeps `[O-]n1cccc1` an anion while repairing `[O-]n1ccccc1`.

def test_pyridine_n_oxide_written_with_a_double_bond():
    # `O=n1ccccc1`.  The neutral N cannot hold that double bond, so it becomes `[O-][n+]1ccccc1` and
    # the ring is six must-match atoms.  Left as written it is five carbons around a saturated N --
    # an odd count with nowhere to go
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'O'], cycle(6), extra=[(5, 6, 2)])
    assert kekule_classify(mol, bonds)[ids[5]] == 'must'
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and unresolved == []
    assert len(log) == 1 and 'cannot carry a double bond' in log[0]
    assert (mol.charge_of(ids[5]), mol.charge_of(ids[6])) == (1, -1)
    assert mol.order_of(ids[5], ids[6]) == 1
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_pyridine_n_imide_written_with_a_double_bond():
    # `CN=n1ccccc1` -- the same repair where the exocyclic N carries a substituent, so a one- OR
    # two-coordinate N is the shape and the methyl is part of the case
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'N', 'C'], cycle(6), extra=[(5, 6, 2), (6, 7, 1)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and unresolved == []
    assert len(log) == 1 and 'cannot carry a double bond' in log[0]
    assert (mol.charge_of(ids[5]), mol.charge_of(ids[6])) == (1, -1)
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_the_unconditional_repair_conserves_total_charge():
    # the property that makes it safe to apply with no context: one bond order unit becomes one
    # charge unit on each atom, so neither the molecule's charge nor either atom's valence moves
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'O'], cycle(6), extra=[(5, 6, 2)])
    before = sum(mol.charge_of(i) for i in ids)
    kekule(mol, bonds)
    assert sum(mol.charge_of(i) for i in ids) == before == 0


def test_the_correctly_spelled_n_oxide_is_not_touched():
    # `[O-][n+]1ccccc1` is the form everything above normalises TO, so it must arrive as a no-repair
    # case -- an empty log, not a log line saying the file agreed with itself
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'O'], cycle(6), extra=[(5, 6, 1)],
                               charges={5: 1, 6: -1})
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_pyrazine_dioxide_mis_spelled_in_both_directions():
    # `[O-]n1cc[n+](=O)cc1` -- N1 written neutral with an anion, N4 written cationic with a neutral
    # oxide.  Each half alone would move the total charge; together they cancel, so this is repaired
    # without needing the ring to fail first.  It WOULD otherwise kekulise: both nitrogens saturated
    # leaves four carbons in two adjacent pairs, a valence-legal 1,4-dihydropyrazine that throws
    # away the aromaticity the input asserted on all six bonds
    mol, ids, bonds, _ = build(['N', 'C', 'C', 'N', 'C', 'C', 'O', 'O'], cycle(6),
                               extra=[(0, 6, 1), (3, 7, 2)], charges={3: 1, 6: -1})
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and unresolved == []
    assert len(log) == 2 and all('mis-spelled in both directions' in line for line in log)
    assert [mol.charge_of(i) for i in ids] == [1, 0, 0, 1, 0, 0, -1, -1]
    assert sum(mol.charge_of(i) for i in ids) == 0
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_furoxan_written_with_a_neutral_n_oxide():
    # `O=[n+]1ccno1`, 1,2,5-oxadiazole 2-oxide.  One mis-spelling, no partner to cancel it, so the
    # repair costs the molecule its spurious +1 -- and is taken only because the ring cannot
    # kekulise as written.  Ring O donates and the two carbons pair off, so the oxide nitrogen and
    # the plain ring N are two must-match atoms with nothing to match.
    #
    # STATED H0 ON THE PLAIN RING NITROGEN IS PART OF THE CASE, not fixture noise: without it that N
    # is the pyrrole/pyridine free choice, it absorbs the odd count as a donor, and the system
    # kekulises with the mis-spelled oxide left standing.  Bare `n` in SMILES states zero, which is
    # why the string fails where a builder that says nothing does not
    mol, ids, bonds, h = build(['N', 'C', 'C', 'N', 'O', 'O'], cycle(5), extra=[(0, 5, 2)],
                               charges={0: 1}, stated_h={3: 0})
    changed, log, unresolved = spread(kekule(mol, bonds, h))
    assert changed and unresolved == []
    assert len(log) == 1 and 'has no Kekule form as written' in log[0]
    assert (mol.charge_of(ids[0]), mol.charge_of(ids[5])) == (1, -1)
    assert mol.order_of(ids[0], ids[5]) == 1
    assert sum(mol.charge_of(i) for i in ids) == 0


def test_pyridine_n_olate_in_a_six_ring_gains_the_cation_it_needs():
    # `[O-]n1ccccc1`: a neutral three-coordinate N donating into a six-ring leaves five carbons, so
    # there is no Kekule form and the only reading that gives one is pyridine N-oxide
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'O'], cycle(6), extra=[(5, 6, 1)],
                               charges={6: -1})
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and unresolved == []
    assert len(log) == 1 and 'has no Kekule form as written' in log[0]
    assert (mol.charge_of(ids[5]), mol.charge_of(ids[6])) == (1, -1)
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_the_same_shape_in_a_five_ring_keeps_its_charge():
    # `[O-]n1cccc1` is pyrrol-1-olate: the same neutral N-with-an-anion, and here it is simply
    # right -- a pyrrole-type nitrogen has three sigma bonds and the four carbons pair off.  The
    # molecule stays an anion, and no rule in this file mentions ring size to make that happen
    mol, ids, bonds, _ = build(['C'] * 4 + ['N', 'O'], cycle(5), extra=[(4, 5, 1)], charges={5: -1})
    assert clean(kekule(mol, bonds))
    assert mol.charge_of(ids[4]) == 0
    assert sum(mol.charge_of(i) for i in ids) == -1


def test_the_n_oxide_ylide_spelling_is_left_alone():
    # `O=[n+]1cccc[c-]1` -- a cationic N holding a double bond, balanced by a ring carbanion.  Every
    # atom is valence-legal, the total charge is zero and the ring HAS a Kekule form, so the
    # unbalanced gate never opens.  Repairing it would hand back an anion
    mol, ids, bonds, _ = build(['N'] + ['C'] * 5 + ['O'], cycle(6), extra=[(0, 6, 2)],
                               charges={0: 1, 5: -1})
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and unresolved == [] and log == []
    assert (mol.charge_of(ids[0]), mol.charge_of(ids[6])) == (1, 0)
    assert mol.order_of(ids[0], ids[6]) == 2


def test_thiopyranium_olate_has_no_kekule_form_and_says_so():
    # `[O-][s+]1ccccc1`.  THE ONE DELIBERATE REFUSAL in this section: a three-coordinate cationic S
    # takes no ring double bond, which leaves five carbons -- an odd count, so no Kekule form exists
    # for any spelling of the substituent.  chython 2 has a rule for this shape that rewrites it to
    # `S=O`; that is not ported, because the neutral three-coordinate S it destroys is the only
    # state that could ever have matched, and the rewrite loses the input's charges as well
    mol, ids, bonds, _ = build(['C'] * 5 + ['S', 'O'], cycle(6), extra=[(5, 6, 1)],
                               charges={5: 1, 6: -1})
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == [tuple(sorted(ids[:6]))]
    assert len(log) == 1 and 'no Kekule form' in log[0]
    # the refusal costs nothing: the charges are the ones the input stated, and the ring came back
    # with the best partial assignment rather than an exception
    assert (mol.charge_of(ids[5]), mol.charge_of(ids[6])) == (1, -1)


def test_a_phosphinine_oxide_keeps_its_double_bond():
    # `O=p1ccccc1`.  The repair is nitrogen-only, and this is why: P really can be pentavalent, so
    # the classifier's P arm already accepts the exocyclic double bond and the ring is five carbons
    # around a saturated P.  That has no Kekule form either, but rewriting the phosphorus would be
    # asserting a charge separation that the chemistry does not require
    mol, ids, bonds, _ = build(['C'] * 5 + ['P', 'O'], cycle(6), extra=[(5, 6, 2)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == [tuple(sorted(ids[:6]))]
    assert mol.charge_of(ids[5]) == 0 and mol.order_of(ids[5], ids[6]) == 2


def test_a_rejected_balanced_pair_leaves_the_matching_it_found_intact():
    """A balanced pair is TRIED on a system that already kekulised, so a failed try must cost nothing.

    A benzene fused to an eight-ring whose two mis-spelled nitrogens are walled in by exocyclic
    carbonyls: the pair is balanced, so it is offered, and promoting either nitrogen to a must-match
    atom leaves it with no partner that could take a double bond.  The trial fails, and the benzene's
    three double bonds -- found before the trial reset the matching to look for a better one -- have
    to come back exactly as they were.  Nothing is logged, because nothing happened.
    """
    elements = ['C'] * 7 + ['N', 'C', 'N', 'C', 'C'] + ['O'] * 6
    #           0-5 benzene, 6 C=O, 7 N with an anion, 8 C=O, 9 cationic N with a neutral oxide,
    #           10 C=O, 11 C=O; the eight-ring shares the 0-1 bond with the benzene
    ring6 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    ring8 = [(1, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 0)]
    mol, ids, bonds, _ = build(elements, ring6 + ring8,
                               extra=[(7, 12, 1), (9, 13, 2), (6, 14, 2), (8, 15, 2), (10, 16, 2),
                                      (11, 17, 2)],
                               charges={9: 1, 12: -1})
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert changed and log == [] and unresolved == []
    assert [mol.charge_of(i) for i in ids] == [0] * 9 + [1, 0, 0, -1] + [0] * 5
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_the_six_vendor_spellings_end_to_end():
    """The six strings as they arrived, through the reader, in one table.

    THE ONLY PARSED FIXTURES IN THIS FILE, and the exception is earned twice over.  These are
    verbatim vendor input, so a string is the shape they actually have -- and the reader's hydrogen
    convention is load-bearing for two of them: bare `n` states zero hydrogens, which is what makes
    furoxan's ring nitrogen a must-match atom and the system fail.  A builder fixture that says
    nothing about hydrogens kekulises the same molecule without repairing it, which is correct for
    what it was asked and not what the file said.  Every case above pins the mechanism; this pins
    the answer a user gets.
    """
    from chython.core._core import read_smiles

    for string, expect_unresolved, charge in [
        ('O=n1ccccc1', False, 0),             # pyridine N-oxide, double bond on a neutral N
        ('CN=n1ccccc1', False, 0),            # its N-imide analogue
        ('[O-][s+]1ccccc1', True, 0),         # no Kekule form for any spelling -- see the S test
        ('[O-]n1cc[n+](=O)cc1', False, 0),    # pyrazine 1,4-dioxide, wrong in both directions
        ('O=[n+]1cccc[c-]1', False, 0),       # a valid ylide; left exactly as written
        ('O=[n+]1ccno1', False, 0),           # furoxan, spelled with a spurious +1
    ]:
        mol = read_smiles(string)
        result = mol.kekule()
        ids = [n for n in mol]
        assert bool(result.unresolved) is expect_unresolved, (string, result.log)
        assert sum(mol.charge_of(n) for n in ids) == charge, (string, result.log)
        # whatever it decided, it decided once: the second call finds nothing left to do
        assert quiet(mol.kekule()), (string, result.log)


# --- the hydrogen relaxation: an aromatic nitrogen whose count the notation never carried
#
# PARSED FIXTURES, and for the reason the vendor table above gives: the whole defect is that a bare
# lowercase `n` reaches the arena with a manufactured zero, so the ring is a must-match atom short
# of a Kekule form.  A builder fixture that says nothing about hydrogens stores H_UNKNOWN, which is
# the free choice, and kekulises these rings without ever reaching the repair -- correct for what it
# was asked, and not the question.  `test_stated_nh0_in_a_five_ring_is_repaired_with_the_hydrogen_
# it_needs` is the builder half, with the zero stated outright.

def test_the_opensmiles_five_rings_are_repaired_with_the_hydrogen_the_string_omitted():
    """The common case in the wild, all four of them.

    OpenSMILES is explicit that a lowercase ring `n` decides its hydrogen count from the ring: with
    one it donates a lone pair and takes no double bond, without one it contributes a single pi
    electron and must take one.  Toolkits write the two-electron form without the hydrogen anyway --
    `c1cncn1` for imidazole is what an aromatic SMILES of imidazole usually looks like -- and none
    of these strings has a Kekule form as written.  Refusing them would mean refusing a large share
    of the aromatic SMILES in circulation.

    Compared against the SAME MOLECULE spelled with the hydrogen, not against a fixed string: which
    of two equivalent tautomers the search reaches is a free choice, and for these rings the two are
    the same molecule anyway.
    """
    from chython.core._core import read_smiles

    for bare, spelled in [
        ('c1ccnc1', 'c1cc[nH]c1'),            # pyrrole
        ('c1cncn1', 'c1cnc[nH]1'),            # imidazole
        ('c1cnnc1', 'c1cn[nH]c1'),            # pyrazole
        ('c1nnnn1', 'c1nn[nH]n1'),            # 1H-tetrazole
    ]:
        mol = read_smiles(bare)
        result = mol.kekule()
        assert result.unresolved == [], (bare, result.log)
        assert any('one hydrogen' in line for line in result.log), (bare, result.log)
        assert sum(mol.implicit_h_of(n) or 0 for n in mol if mol.atom(n).element == 7) == 1, bare
        reference = read_smiles(spelled)
        reference.kekule()
        assert mol == reference, (bare, mol.smiles, reference.smiles)
        # and it decided once: the stored orders and counts are already the answer
        assert quiet(mol.kekule()), (bare, result.log)


def test_a_fused_ring_gets_its_hydrogen_on_the_nitrogen_that_needs_it():
    """Benzimidazole and purine, each written with neither of its two nitrogens carrying a hydrogen.

    ONE hydrogen, not two, and that is the property worth testing: the relaxation hands every
    candidate to the search at once and does not enumerate subsets, so the count comes out of the
    matching rather than out of a preference.  Both nitrogens go free, the search uses whichever one
    the ring forces it to use, and only the leftover is protonated.

    WHICH one it is stays a free choice, so both tautomers are accepted.  For benzimidazole the two
    are the same molecule by symmetry; for purine they are not, and picking one here would pin a
    search order rather than a chemical fact.  The input did not say which tautomer it meant, and
    neither does the answer.
    """
    from chython.core._core import read_smiles

    for bare, spellings in [('c1ccc2c(c1)ncn2', ['c1ccc2[nH]cnc2c1', 'c1ccc2nc[nH]c2c1']),
                            ('c1cnc2ncnc2c1', ['c1cnc2[nH]cnc2c1', 'c1cnc2nc[nH]c2c1'])]:
        mol = read_smiles(bare)
        result = mol.kekule()
        assert result.unresolved == [], (bare, result.log)
        assert sum('one hydrogen' in line for line in result.log) == 1, (bare, result.log)
        references = []
        for spelled in spellings:
            reference = read_smiles(spelled)
            reference.kekule()
            references.append(reference)
        assert any(mol == reference for reference in references), \
            (bare, mol.smiles, [r.smiles for r in references])


def test_a_ring_that_kekulises_as_written_is_never_handed_a_hydrogen():
    """The gate that makes the repair free: pyridine and its relatives never reach it.

    Each of these has a Kekule form with the counts exactly as the reader stored them, so the first
    search succeeds and no relaxation is offered.  Without the gate the nitrogens here are the same
    candidates as the ones above, and pyridine would come back as a 1,2-dihydropyridine radical.
    """
    from chython.core._core import read_smiles

    for string in ['c1ccncc1',                # pyridine
                   'c1cncnc1',                # pyrimidine
                   'c1ccnnc1',                # pyridazine
                   'c1cc[nH+]cc1',            # pyridinium, the hydrogen already stated
                   'c1ccc2ncccc2c1',          # quinoline
                   'c1cc[n-]c1',              # pyrrolide, an anion and not a candidate
                   'c1ccoc1',                 # furan: O is unambiguous, so there is no choice
                   'c1ccsc1']:                # thiophene, likewise
        mol = read_smiles(string)
        result = mol.kekule()
        assert result.unresolved == [], (string, result.log)
        assert not any('one hydrogen' in line for line in result.log), (string, result.log)
        for n in mol:
            if mol.atom(n).element == 7:
                assert mol.implicit_h_of(n) == read_smiles(string).implicit_h_of(n), (string, n)


# A STATED HYDROGEN IS NOT PRIVILEGED OVER A STATED ZERO.  `c1cc[nH]cc1` does not come back
# unresolved: the surplus-hydrogen relaxation gives the hydrogen up instead -- see
# `test_a_stated_nh_in_a_six_ring_gives_the_hydrogen_up_rather_than_being_refused` at the end of this
# file.


def test_repairing_an_n_oxide_twice_is_the_same_as_repairing_it_once():
    # the fixed point of every rule here is the charge-separated form, so a second call has nothing
    # to find: no charge moves and the log is empty
    mol, ids, bonds, _ = build(['C'] * 5 + ['N', 'O'], cycle(6), extra=[(5, 6, 2)])
    kekule(mol, bonds)
    charges = [mol.charge_of(i) for i in ids]
    orders_before = orders(mol, bonds)
    assert quiet(kekule(mol, bonds))
    assert [mol.charge_of(i) for i in ids] == charges
    assert orders(mol, bonds) == orders_before


# --- the repairs: each is a log line and a molecule, never an exception

def test_aromatic_bond_in_no_ring_is_read_single():
    # MDL fixture 1: bond type 4 on a bond that is in no ring
    mol, ids, bonds, _ = build('C' * 7, cycle(6) + [(0, 6)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == []
    assert len(log) == 1 and 'in no ring' in log[0]
    assert mol.order_of(ids[0], ids[6]) == 1
    assert alternating(mol, bonds[:6]) == {i: 1 for i in ids[:6]}


def test_biphenyl_written_as_one_aromatic_system():
    # `c1ccccc1c2ccccc2` -- the inter-ring bond is on no cycle OF THE MOLECULE, so it is single
    # however it was written, and both rings still kekulise.  This is the case that stops the
    # ring-membership test below from over-reaching: an aromatic bond between two rings looks
    # exactly like an aromatic bond inside one until you ask the right graph.
    aromatic = cycle(6) + cycle(6, 6) + [(0, 6)]
    mol, ids, bonds, _ = build('C' * 12, aromatic)
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == []
    assert len(log) == 1 and 'in no ring' in log[0]
    assert mol.order_of(ids[0], ids[6]) == 1
    assert alternating(mol, bonds) == {i: 1 for i in ids}


# --- an explicit single bond inside an aromatic ring is a preference, not a wall
#
# Ring membership is a property of the MOLECULE.  An aromatic bond that is a bridge of the stated
# aromatic SUBGRAPH may still sit inside a ring, and then it carries ring pi electrons like any
# other.  Every case below is a real string chython 2 accepts, and chython 2's answer to all three
# is benzene -- named here as the oracle because this trio was a live wrong answer: the subgraph
# bridge test read every one of them as cyclohexane.

def test_one_explicit_single_bond_in_an_aromatic_ring():
    # `c1c-cccc1` -- five stated aromatic bonds forming a path, every one a bridge of that path
    # and every one inside the same six-ring.  chython 2: C1=CC=CC=C1
    aromatic = [(0, 1), (2, 3), (3, 4), (4, 5), (5, 0)]
    mol, ids, bonds, _ = build('C' * 6, aromatic, extra=[(1, 2, 1)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}
    assert mol.order_of(ids[1], ids[2]) == 1        # the stated single bond is respected


def test_two_explicit_single_bonds_in_an_aromatic_ring():
    # `c1cc-cc-c1` -- four stated aromatic bonds, a three-atom path plus a lone edge.
    # chython 2: C1=CC=CC=C1
    aromatic = [(0, 1), (1, 2), (3, 4), (5, 0)]
    mol, ids, bonds, _ = build('C' * 6, aromatic, extra=[(2, 3, 1), (4, 5, 1)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_alternating_explicit_single_bonds_in_an_aromatic_ring():
    # `c1c-cc-cc-1` -- three stated aromatic bonds, pairwise disjoint, so each one is forced
    # double and the answer is benzene without any search.  chython 2: C1=CC=CC=C1
    aromatic = [(0, 1), (2, 3), (4, 5)]
    mol, ids, bonds, _ = build('C' * 6, aromatic, extra=[(1, 2, 1), (3, 4, 1), (5, 0, 1)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_partially_aromatic_ring():
    # MDL fixture 2: a vendor file marks four bonds of a six-ring aromatic and leaves two alone,
    # so the marked run spans FIVE atoms -- `c1cccc-c-1`.  Those five are all must-match carbons
    # and five cannot pair off, so the honest answer is "no Kekule form" with a maximal partial
    # assignment, not silence.  chython 2 refuses this one outright, with
    # `InvalidAromaticRing: not in ring aromatic bond or hypercondensed rings: {1, 2, 3, 4, 5}` --
    # it raises whenever a ring atom carries no aromatic bond at all, because its promotion step
    # (see `arom_prune_acyclic`) demands that every atom of the ring be aromatic-bonded.  Both
    # implementations reject; only one of them does it without an exception.
    aromatic = [(0, 1), (1, 2), (2, 3), (3, 4)]
    mol, ids, bonds, _ = build('C' * 6, aromatic, extra=[(4, 5, 1), (5, 0, 1)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == [tuple(ids[:5])]
    assert len(log) == 1 and 'no Kekule form' in log[0]
    assert 'leaves 1 atom(s) with an incomplete valence' in log[0]
    assert mol.order_of(ids[4], ids[5]) == 1
    assert mol.order_of(ids[5], ids[0]) == 1


def test_partially_aromatic_fused_ring():
    # the same idea where the aromatic part IS a ring: one ring of a fused pair is marked
    # aromatic and the other is drawn Kekule.  The marked ring must still alternate
    aromatic = cycle(6)
    drawn = [(4, 6, 1), (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 3, 1)]
    mol, ids, bonds, _ = build('C' * 10, aromatic, extra=drawn)
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}


def test_odd_carbon_ring_has_no_kekule_form():
    # neutral C5: five must-match atoms cannot pair off, and no search order changes that.  The
    # answer is "impossible", reported once, with a maximal partial assignment left behind
    mol, ids, bonds, _ = build('C' * 5, cycle(5))
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert unresolved == [tuple(ids)]
    assert len(log) == 1 and 'no Kekule form' in log[0]
    assert 'leaves 1 atom(s) with an incomplete valence' in log[0]
    counts = alternating(mol, bonds)
    assert sorted(counts.values()) == [1, 1, 1, 1]      # four of five saturated: maximal


def test_stated_nh0_in_a_five_ring_is_repaired_with_the_hydrogen_it_needs():
    """`[nH0]1cccc1`: the input insists the N takes a ring double bond, and four carbons cannot
    let it.  The classification still says so -- and that much chython 2 cannot even represent,
    because its sentinel collapses `[nH0]` and `n` -- but the answer is pyrrole.

    A STATED ZERO IS NOT PRIVILEGED.  The reading that a stated zero is a fact the kekuliser has no
    business overruling does not survive the question "what molecule is `[nH0]1cccc1`?" --
    there is none, a two-coordinate neutral nitrogen with no hydrogen and no double bond is not a
    valence state -- so honouring the zero buys a half-assigned ring and a hypovalent atom instead
    of the one molecule the input could have meant.  A stated zero on a system with no Kekule form
    is garbage input like any other, and repairing it at this boundary is what the rest of the file
    already does with a mis-drawn N-oxide.  A stated hydrogen is not privileged either, and for the
    same reason: see `test_a_stated_nh_in_a_six_ring_gives_the_hydrogen_up_rather_than_being_refused`.
    """
    mol, ids, bonds, h = build(['C'] * 4 + ['N'], cycle(5), stated_h={4: 0})
    assert kekule_classify(mol, bonds, h)[ids[4]] == 'must'
    changed, log, unresolved = spread(kekule(mol, bonds, h))
    assert unresolved == []
    assert len(log) == 1 and 'one hydrogen' in log[0] and str(ids[4]) in log[0]
    # the four carbons pair off and the nitrogen carries the hydrogen the ring needs
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}
    assert mol.implicit_h_of(ids[4]) == 1


def test_non_aromatic_element_is_logged():
    mol, ids, bonds, _ = build(['C'] * 5 + ['Si'], cycle(6))
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert any('has no aromatic form' in line and 'Si' in line for line in log)
    # five carbons around a saturated Si cannot pair off, so this is unresolved too -- the point
    # is that both facts are reported and neither raises
    assert unresolved == [tuple(ids)]


def test_invalid_aromatic_state_is_logged_not_raised():
    # an aromatic N with two stated hydrogens: chython 2 raises `InvalidAromaticRing`
    mol, ids, bonds, h = build(['C'] * 4 + ['N'], cycle(5), stated_h={4: 2})
    changed, log, unresolved = spread(kekule(mol, bonds, h))
    assert any('not a valid aromatic state' in line and '2H' in line for line in log)
    assert unresolved == []
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_triple_bond_into_an_aromatic_ring_is_logged():
    mol, ids, bonds, _ = build(['C'] * 8, cycle(6), extra=[(0, 6, 1), (6, 7, 3)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert (log, unresolved) == ([], [])       # a nitrile SUBSTITUENT is perfectly ordinary
    assert alternating(mol, bonds) == {i: 1 for i in ids[:6]}
    # ... but a triple bond ON a ring atom is not
    mol, ids, bonds, _ = build(['C'] * 6 + ['N'], cycle(6), extra=[(0, 6, 3)])
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert any('triple bond into an aromatic ring' in line for line in log)


def test_hypercondensed_atom_is_logged():
    # four aromatic bonds on one atom cannot happen in a real aromatic system.  Two rings
    # sharing a single atom, not a bond -- spiro, so nothing here is a bridge and the degree
    # survives the pruning pass
    aromatic = cycle(4) + [(0, 4), (4, 5), (5, 6), (6, 0)]
    mol, ids, bonds, _ = build('C' * 7, aromatic)
    changed, log, unresolved = spread(kekule(mol, bonds))
    assert any('hypercondensed' in line for line in log)


# --- contract and hygiene

def test_every_system_that_comes_out_says_so():
    """The kekulisation itself is the event `mol.log` exists to hold, so it is a record.

    ONE PER SYSTEM THE PASS SOLVED, `INFO`, naming that system's atoms in sorted stable ids -- two
    disconnected rings are two lines, and biphenyl is two as well, because the bond between them is
    read single first and the system splits.  A system it could not solve contributes none: its own
    `no-kekule-form` line already says what was written there.

    The counted-against property is the second half: nothing moved means nothing said, which is what
    keeps `canonicalize()`'s fixed-point loop from growing a line per round for a molecule at rest.
    """
    mol, ids, bonds, _ = build('C' * 12, cycle(6) + cycle(6, 6))
    result = kekule(mol, bonds)
    assert [(x.rule, x.atoms, x.severity) for x in result.log] == \
        [('kekule:kekulized', tuple(sorted(ids[:6])), 'info'),
         ('kekule:kekulized', tuple(sorted(ids[6:])), 'info')]
    assert quiet(kekule(mol, bonds)), 'a system already in a Kekule form said something anyway'

    # biphenyl: one stated system, two solved ones
    mol, ids, bonds, _ = build('C' * 12, cycle(6) + cycle(6, 6) + [(0, 6)])
    result = kekule(mol, bonds)
    assert [x.atoms for x in result.log if x.rule == 'kekule:kekulized'] == \
        [tuple(sorted(ids[:6])), tuple(sorted(ids[6:]))]

    # an unresolved thiopyranium olate beside a benzene: the ring that came out is the only one named
    mol, ids, bonds, _ = build(['C'] * 5 + ['S', 'O'] + ['C'] * 6, cycle(6) + cycle(6, 7),
                               extra=[(5, 6, 1)], charges={5: 1, 6: -1})
    result = kekule(mol, bonds)
    assert result.unresolved == [tuple(sorted(ids[:6]))]
    assert [x.atoms for x in result.log if x.rule == 'kekule:kekulized'] == [tuple(sorted(ids[7:]))]


def test_radical_aromatic_heteroatom_needs_no_hydrogen_count():
    """A caller may kekulise BEFORE it has worked out any hydrogen count, and `kekule` writes none.

    The classifier DOES read the stored nibble -- `arom_setup` takes the store as its default source --
    and the affordance holds anyway: a molecule whose counts are not worked out yet holds `H_UNKNOWN`,
    the classifier recognises that as unstated, and answers exactly as it does for a count forced
    unstated.

    The second assertion is the load-bearing one and it is a NEGATIVE: kekulisation must not fill a
    count in as a side effect.  It reads `is None` and not `== 0` because "nobody counted" is the
    sentinel, and `== 0` cannot tell a preserved zero from a written one.
    """
    mol, ids, bonds, _ = build(['C'] * 4 + ['N'], cycle(5), radicals=(4,))
    assert kekule_classify(mol, bonds)[ids[4]] == 'must_not'
    assert all(mol.implicit_h_of(i) is None for i in ids)
    assert clean(kekule(mol, bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids[:4]}


def test_empty_aromatic_set_is_a_no_op():
    mol, ids, bonds, _ = build('CC', [], extra=[(0, 1, 1)])
    assert quiet(kekule(mol, [], None))
    assert mol.order_of(ids[0], ids[1]) == 1


def test_duplicate_and_reversed_pairs_are_accepted():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    doubled = bonds + [(b, a) for a, b in bonds]
    changed, log, unresolved = spread(kekule(mol, doubled))
    assert (log, unresolved) == ([], [])
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_running_twice_is_the_same_as_running_once():
    mol, ids, bonds, _ = build(['C'] * 5 + ['N'], cycle(6))
    assert kekule(mol, bonds).changed
    first = orders(mol, bonds)
    changed = kekule(mol, bonds).changed
    # the second call still assigns, and still assigns the same thing -- but reports that nothing
    # moved, which is the difference between idempotent and merely stable
    assert not changed
    assert orders(mol, bonds) == first


def test_the_method_and_the_function_are_the_same_operation():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    assert clean(mol.kekule(bonds))
    assert alternating(mol, bonds) == {i: 1 for i in ids}


def test_no_argument_on_a_molecule_with_no_aromatic_bonds_is_a_no_op_that_says_so():
    # `build` writes order 1 on the "aromatic" edges, so this molecule stores nothing aromatic and
    # the derived set is empty.  Empty means nothing moved AND nothing to report -- not an error
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    assert mol.aromatic_bond_count == 0
    assert quiet(mol.kekule())
    assert orders(mol, bonds) == [1] * 6


def test_the_no_argument_form_kekulises_the_stored_aromatic_bonds():
    """The ordinary call, now that order 4 is a stored order (it replaces the spent gate tripwire).

    A caller who wants a Kekule form does not have to tell the molecule which of its own bonds are
    aromatic. `arom_stored_bonds` reads them off the half-edges, so this is the same operation as
    passing the set explicitly -- asserted here against the explicit form on the same shape rather
    than against a literal, so the two spellings cannot drift apart.
    """
    mol = MoleculeContainer()
    ids = [mol.add_atom('C') for _ in range(6)]
    pairs = [(ids[i], ids[(i + 1) % 6]) for i in range(6)]
    for a, b in pairs:
        mol.add_bond(a, b, 4)
    assert mol.aromatic_bond_count == 6
    assert clean(mol.kekule())
    assert orders(mol, pairs) == [2, 1, 2, 1, 2, 1]
    assert mol.aromatic_bond_count == 0, 'a kekulised bond is no longer aromatic, flag included'
    assert quiet(mol.kekule()), 'and the second call has nothing left to do'

    explicit = MoleculeContainer()
    eids = [explicit.add_atom('C') for _ in range(6)]
    epairs = [(eids[i], eids[(i + 1) % 6]) for i in range(6)]
    for a, b in epairs:
        explicit.add_bond(a, b, 4)
    assert clean(explicit.kekule(epairs))
    assert orders(explicit, epairs) == orders(mol, pairs)


def test_a_stated_subset_leaves_the_rest_of_the_stored_aromatic_bonds_alone():
    """Why the explicit form survives the arena learning to store order 4.

    Biphenyl written aromatic, with only the first ring stated: the second ring keeps order 4 and the
    count says so. A caller that could not do this would have no way to kekulise a vendor file's
    per-bond marks without touching bonds the file did not mark.
    """
    mol = MoleculeContainer()
    ids = [mol.add_atom('C') for _ in range(12)]
    first = [(ids[i], ids[(i + 1) % 6]) for i in range(6)]
    second = [(ids[6 + i], ids[6 + (i + 1) % 6]) for i in range(6)]
    for a, b in first + second:
        mol.add_bond(a, b, 4)
    mol.add_bond(ids[0], ids[6], 1)
    assert clean(mol.kekule(first))
    assert mol.aromatic_bond_count == 6, 'the unstated ring is untouched'
    assert orders(mol, second) == [4] * 6


def test_a_bond_that_is_not_a_bond_raises_key_error():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    with raises(KeyError):
        kekule(mol, [(ids[0], ids[3])])
    with raises(KeyError):
        kekule(mol, [(ids[0], 999)])


def test_a_self_loop_raises_value_error():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    with raises(ValueError):
        kekule(mol, [(ids[0], ids[0])])


def test_pending_edits_are_refused():
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    with raises(RuntimeError):
        with mol.edit():
            mol.add_atom('C')
            kekule(mol, bonds)


def test_the_result_is_a_named_object_and_not_a_tuple():
    # this return value already grew from two fields to three.  A caller that unpacks it
    # positionally is the call site that breaks when it grows again, so unpacking must not work at
    # all -- a TypeError today is cheaper than a silent misread later
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    result = mol.kekule(bonds)
    assert isinstance(result, KekuleResult)
    assert not isinstance(result, tuple)
    with raises(TypeError):
        a, b, c = result
    assert 'changed=True' in repr(result)


def test_kekule_copy_leaves_its_input_alone():
    # the shape a format bridge needs: it is holding somebody else's molecule, and input fidelity
    # is an invariant of `molecule_to_inchi`.  Registered with the InChI side at module init
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    out = kekule_copy(mol, bonds)
    assert orders(mol, bonds) == [1] * 6
    assert alternating(out, bonds) == {i: 1 for i in ids}


def test_kekule_copy_raises_rather_than_hand_over_wrong_orders():
    # a partial assignment is the right answer for a caller who can read the log, and the wrong
    # answer for libinchi, which would return a confidently wrong InChI for it
    mol, ids, bonds, _ = build('C' * 5, cycle(5))
    with raises(AromaticKekulizeError, match='no Kekule form'):
        kekule_copy(mol, bonds)
    assert orders(mol, bonds) == [1] * 5


def test_the_error_type_exists_and_is_distinguishable():
    # a caller stubbing against this entry point needs to tell a real bug from an input it should
    # have expected, which is the whole reason this type is not ValueError
    assert issubclass(AromaticKekulizeError, RuntimeError)
    assert not issubclass(AromaticKekulizeError, ValueError)


def test_the_harness_can_fail():
    # ruling F102: the clean sweep above is evidence only if `alternating` could have failed.
    # Kekulise benzene, then break one bond back to single by hand -- the check must notice
    mol, ids, bonds, _ = build('C' * 6, cycle(6))
    kekule(mol, bonds)
    assert alternating(mol, bonds) == {i: 1 for i in ids}
    for a, b in bonds:
        if mol.order_of(a, b) == 2:
            mol.set_order(a, b, 1)
            break
    assert alternating(mol, bonds) != {i: 1 for i in ids}


# --- the fourth and fifth relaxations: a missing formal charge and a surplus hydrogen


def test_an_n_alkyl_pyridinium_written_neutral_gets_the_charge_it_needs():
    """`c1ccn(C)cc1` -- five carbons and a three-coordinate neutral nitrogen cannot pair off.

    The nitrogen has no room for a ring double bond either way, so the ring is not asking for a
    hydrogen (there is nowhere to put one) and not asking for a bond order: it is asking for the
    cation the substituent already implies.  N-methylpyridinium and N-methylnicotinamide arrive
    from real files spelled this way constantly.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1ccn(C)cc1')
    changed, log, unresolved = spread(mol.kekule())
    assert unresolved == []
    assert changed
    assert any('read as' in line and '(+)' in line for line in log), log
    n = next(i for i in mol if mol.atom(i).element == 7)
    assert mol.charge_of(n) == 1


def test_an_n_hydroxy_pyridine_is_read_as_the_charge_separated_n_oxide():
    """`c1ccn(O)cc1` -- the same failed ring, but with a proton to spare on the substituent.

    Two trials would fix the ring: `[n+](O)` alone, which hands back a cation, and `[n+][O-]` with
    the hydroxyl's proton gone, which is pyridine N-oxide and conserves the molecule's charge.  The
    charge-conserving one is tried first, so that is the one that lands.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1ccn(O)cc1')
    changed, log, unresolved = spread(mol.kekule())
    assert unresolved == []
    n = next(i for i in mol if mol.atom(i).element == 7)
    o = next(i for i in mol if mol.atom(i).element == 8)
    assert mol.charge_of(n) == 1
    assert mol.charge_of(o) == -1
    assert mol.implicit_h_of(o) == 0


def test_an_azole_carrying_a_hydrogen_it_cannot_afford_gives_it_up():
    """`Cn1cc[nH]c1` -- both nitrogens want to be two-electron donors and only one may be.

    A five-ring with two donors has three carbons left over, which is odd, so nothing pairs off.
    The N-methyl nitrogen cannot give anything up; the `[nH]` can, and once it does the ring is
    1-methylimidazole.  This is the case the corpus produced most often after the missing charge.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('Cn1cc[nH]c1')
    changed, log, unresolved = spread(mol.kekule())
    assert unresolved == []
    assert changed
    assert any('surplus' in line for line in log), log
    assert sorted(mol.implicit_h_of(i) for i in mol if mol.atom(i).element == 7) == [0, 0]


def test_a_stated_nh_in_a_six_ring_gives_the_hydrogen_up_rather_than_being_refused():
    """`c1cc[nH]cc1`, which chython 2 answers with `InvalidAromaticRing`.

    A STATED HYDROGEN IS NOT PRIVILEGED EITHER.  The reading that a hydrogen the input DID state is a
    fact the kekuliser has no business taking away does not survive the question the file already
    answers for a stated zero: what molecule is `c1cc[nH]cc1`?  There is none -- a six-ring cannot host
    a two-electron donor and pair off five carbons -- so honouring the hydrogen buys a half-assigned
    ring instead of the one molecule the input could have meant, which is pyridine.

    A stated hydrogen and a stated zero are therefore treated alike, and the symmetry is the
    point: both are garbage input the ring contradicts, and both are repaired at this boundary.  The
    gate is unchanged -- a ring that kekulises as written is never offered either relaxation -- so
    `c1cc[nH]c1` is still pyrrole and keeps its hydrogen.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1cc[nH]cc1')
    changed, log, unresolved = spread(mol.kekule())
    assert unresolved == []
    assert changed
    assert any('surplus' in line for line in log), log
    assert [mol.implicit_h_of(n) for n in mol if mol.atom(n).element == 7] == [0]


def test_neither_new_relaxation_touches_a_ring_that_kekulises_as_written():
    """The gate, and it is the same gate the other three relaxations pass through."""
    from chython.core._core import read_smiles

    for string in ['c1ccncc1',                # pyridine
                   'c1cc[nH]c1',              # pyrrole: the hydrogen is needed, so it stays
                   'c1cc[n-]c1',              # pyrrolide
                   'c1ccoc1',                 # furan
                   'c1ccsc1',                 # thiophene
                   'Cn1cccc1',                # 1-methylpyrrole: the N is already a donor
                   'c1cc[nH+]cc1',            # pyridinium, charge already stated
                   'C[n+]1ccccc1',            # N-methylpyridinium spelled correctly
                   'c1ccc2ncccc2c1',          # quinoline
                   'c1ccc2[nH]ccc2c1']:       # indole
        # NOT `O=[n+]1ccccc1[O-]` here: `arom_separate_charges` rewrites its `[n+]=O` to `[n+][O-]`
        # unconditionally, before any relaxation is offered, so a charge does move on it and the
        # comparison below would fail for a reason that has nothing to do with this test.
        mol = read_smiles(string)
        before = [(mol.charge_of(n), mol.implicit_h_of(n)) for n in mol]
        result = mol.kekule()
        assert result.unresolved == [], (string, result.log)
        assert not any('surplus' in line or '(+)' in line for line in result.log), (string,
                                                                                    result.log)
        assert [(mol.charge_of(n), mol.implicit_h_of(n)) for n in mol] == before, string


def test_a_ring_no_relaxation_can_reach_is_still_reported():
    """Five aromatic carbons: an odd count of must-match atoms, and no candidate of any kind.

    There is no charge and no hydrogen that makes an odd number even, so this is the fixture for
    "reported honestly" now that the azoles are repaired.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1cccc1')
    changed, log, unresolved = spread(mol.kekule())
    assert len(unresolved) == 1 and len(unresolved[0]) == 5
    assert any('no Kekule form' in line for line in log), log


def test_the_new_relaxations_are_idempotent():
    from chython.core._core import read_smiles

    for string in ['c1ccn(C)cc1', 'c1ccn(O)cc1', 'Cn1cc[nH]c1', 'c1cc[nH]cc1']:
        mol = read_smiles(string)
        assert mol.kekule().changed, string
        charges = [mol.charge_of(n) for n in mol]
        hydrogens = [mol.implicit_h_of(n) for n in mol]
        result = mol.kekule()
        assert not result.changed, (string, result.log)
        assert [mol.charge_of(n) for n in mol] == charges, string
        assert [mol.implicit_h_of(n) for n in mol] == hydrogens, string


def test_a_neutral_ring_atom_with_no_hydrogen_is_not_offered_a_charge_it_cannot_carry():
    """The predicate is the classification table, so an arm that answers `may` is not a candidate.

    `c1cc[n]cc1` with a two-coordinate nitrogen at charge +1 classifies MAY, not MUST -- a
    pyridinium nitrogen may take a ring double bond or not -- so the cation relaxation declines it
    and the hydrogen relaxation, which is what that ring actually needs, gets it.  Asserting the
    hydrogen line rather than the charge line is how that ordering is pinned.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1cncn1')            # imidazole with no hydrogen stated anywhere
    changed, log, unresolved = spread(mol.kekule())
    assert unresolved == []
    assert any('one hydrogen' in line for line in log), log
    assert not any('(+)' in line for line in log), log
    assert sum(mol.charge_of(n) for n in mol) == 0


# --- the aromatic bond that is in no ring, and what "read as single" has to mean

def test_biphenyl_written_the_way_opensmiles_requires_it_be_read():
    """`c1ccccc1c1ccccc1` is biphenyl, and the reader is right to make that bond aromatic.

    OpenSMILES is explicit: an unspecified bond between two aromatic atoms is an AROMATIC bond, which
    is why biphenyl has to be written `c1ccccc1-c1ccccc1` to mean a single bond between the rings.  So
    the reader storing order 4 there is the spec being followed, not a defect, and dealing with it is
    this file's job -- the inter-ring bond lies on no cycle, so it is pruned from the aromatic set and
    logged "read as single", and then each ring is an ordinary benzene.

    THE CLASSIFIER MUST NOT READ THE PRUNED BOND AS A DOUBLE.  An `order >= 2` test passes on order 4,
    which makes each ipso carbon look like a quinone carbonyl carbon and answer must-not, leaving five
    must-match atoms in a six-ring; an odd count has no perfect matching, so BOTH rings of a plain
    biphenyl come back unresolved with an atom each left unsaturated.  Two spellings of the commonest
    biaryl in chemistry have to agree.
    """
    from chython.core._core import read_smiles

    implied, explicit = read_smiles('c1ccccc1c1ccccc1'), read_smiles('c1ccccc1-c1ccccc1')
    changed, log, unresolved = spread(implied.kekule())
    assert unresolved == []
    assert any('is in no ring; read as single' in line for line in log), log
    assert not explicit.kekule().unresolved
    assert implied.canonical_bytes == explicit.canonical_bytes, (str(implied), str(explicit))


def test_the_biaryl_bond_is_read_as_single_whatever_it_joins():
    """Same shape across heteroaromatics and two bonds deep, plus the branch spelling.

    The branch form matters on its own: it rules out the ring-closure syntax as the cause, which is
    what a reader would suspect first.  Every ring atom must end up with exactly one ring double
    bond -- that is what "resolved" means, and asserting it here rather than trusting `unresolved`
    keeps the two independent.
    """
    from chython.core._core import read_smiles

    for string in ['c1ccc(cc1)c1ccccc1',                # the branch spelling of biphenyl
                   'c1ccccc1c1ccncc1',                  # 2-phenylpyridine
                   'c1ccccc1c1cc[nH]c1',                # 2-phenylpyrrole
                   'c1ccccc1c1ccccc1c1ccccc1',          # o-terphenyl: two bonds to prune
                   'c1ccccc1c1ccc(cc1)c1ccccc1']:       # p-terphenyl
        mol = read_smiles(string)
        changed, log, unresolved = spread(mol.kekule())
        assert unresolved == [], (string, log)
        for n in mol:
            doubles = sum(1 for m in mol.neighbors_of(n) if mol.order_of(n, m) == 2)
            if mol.element_of(n) == 6 and not mol.implicit_h_of(n) is None:
                assert doubles == 1, (string, n, str(mol))


def test_an_unresolved_system_says_it_was_rewritten_rather_than_left_as_drawn():
    """The message carries the whole claim here, because the behaviour is the intended one.

    A system with no Kekule form gets the best matching there is -- the record stays readable and the
    rest of the molecule stays usable, which is what accepting garbage input requires.  So the message
    must not read as though the system came back as drawn, which "N atoms left unsaturated" alone does.
    The aromatic flags are gone either way, and the atoms the matching could not pair carry
    single bonds and an incomplete valence with NO radical flag -- deliberately, because one failed
    matching is no evidence that the input meant a radical.
    """
    from chython.core._core import read_smiles

    mol = read_smiles('c1cccc1')
    changed, log, unresolved = spread(mol.kekule())
    assert len(unresolved) == 1
    line, = [x for x in log if 'no Kekule form' in x]
    assert 'best matching' in line, line
    assert 'incomplete valence' in line, line
    assert 'no radical flag' in line, line
    assert not any(mol.radical_of(n) for n in mol), 'a radical flag would be a claim about the input'
    assert not any(mol.order_of(n, m) == 4 for n in mol for m in mol.neighbors_of(n)), \
        'the system was rewritten, so no aromatic bond may survive in it'
