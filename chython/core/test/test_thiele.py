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
"""`thiele()`: Kekule bond orders become order 4, where a sextet can be spelled.

THE NEGATIVE CASES ARE HALF OF THIS FILE, and they are not symmetrical with each other -- there are
two ways not to be aromatic and the difference is observable:

* NEVER A CANDIDATE.  Cyclobutadiene (wrong size), cyclohexene (two sp3 atoms), fulvene and every
  quinone and pyridone (a ring atom's double bond points out of the ring).  Nothing is refused and
  nothing is logged: no system was proposed, so none was declined.
* PROPOSED AND DECLINED.  The cyclopentadienyl cation, borole, 1H-azepine -- a Kekule match exists,
  every atom is happy, and the electron count is not 4n+2.  These arrive in `.refused` with a line
  in `.log`, because a caller who wrote a five-ring of sp2 carbons and got no aromatic bonds is owed
  the reason.

Molecules are hand-built rather than parsed: hydrogen counts here are load-bearing (an NH is a donor,
an N with no H and two neighbours is not), so stating every count in the fixture is the point rather
than an inconvenience.

WHERE THIS DISAGREES WITH chython 2, measured, all four rows deliberate and all four in
`docs/superpowers/specs/2026-09-02-smiles-write-design.md` section 14.2:

| molecule                   | V2  | V3  | why                                                     |
|----------------------------|-----|-----|---------------------------------------------------------|
| tropylium                  | no  | YES | V2's element filter has no seven-ring carbocation rule   |
| selenophene, tellurophene  | no  | YES | V2's element list stops at S                             |
| borole                     | YES | no  | 4 pi: neutral B is the EMPTY-orbital atom, not a donor    |
| pyridones, p-benzoquinone  | no  | no  | agreement, by different mechanisms                       |

V2 is an oracle, not a contract; each of these was measured against `chython.smiles(...).thiele()`
before it was written down.
"""
from itertools import permutations
from math import factorial
from random import Random

from pytest import mark, raises

from chython.core import H_UNKNOWN, MoleculeContainer
from chython.core._core import ThieleResult, kekule, thiele


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING.  An atom is `(element, implicit_h)` or `(element, implicit_h, charge)`; a bond
# is `(i, j, order)` over atom indices; `order` is a creation order, so every claim can be swept.
def build(atoms, bonds, order=None):
    mol = MoleculeContainer()
    ids = {}
    for k in (range(len(atoms)) if order is None else order):
        spec = atoms[k]
        ids[k] = mol.add_atom(spec[0], implicit_h=spec[1],
                              charge=spec[2] if len(spec) > 2 else 0)
    for i, j, o in bonds:
        mol.add_bond(ids[i], ids[j], o)
    return mol, ids


def ring(n, pattern, offset=0):
    """The bonds of a cycle `offset .. offset + n - 1`, `pattern[i]` the order of bond i -> i + 1."""
    return [(offset + i, offset + (i + 1) % n, 2 if pattern[i] == '=' else 1) for i in range(n)]


def cyc(seq):
    """The undirected edges of a cycle through `seq`, as index pairs."""
    seq = list(seq)
    return {frozenset(p) for p in zip(seq, seq[1:] + seq[:1])}


# The two records a successful pass writes about the representation it rewrote: `thiele()` says which
# ring system it made aromatic, `kekule()` which system it wrote a Kekule form of.  Both are INFO and
# both are the WORK, which is why `mol.log` is not empty after a pass that repaired nothing -- and why
# every assertion below filters them out, since each case here asks what its input NEEDED reporting.
# The records themselves are pinned by `test_every_system_that_comes_out_says_so` and by
# `core/test/test_container_log.py`.
ROUTINE = frozenset(('thiele:aromatized', 'kekule:kekulized'))


def notices(result):
    """The result's log without the routine line for the representation it rewrote."""
    return [x for x in result.log if x.rule not in ROUTINE]


def aromatic(mol, ids, bonds):
    """Which bonds of the fixture now hold order 4, back in fixture indices."""
    return {frozenset((i, j)) for i, j, _ in bonds if mol.order_of(ids[i], ids[j]) == 4}


def refused(mol, ids, result):
    """`.refused`, back in fixture indices, as a set of frozensets."""
    back = {sid: i for i, sid in ids.items()}
    return {frozenset(back[sid] for sid in system) for system in result.refused}


ALL = 'all'      # every bond of the fixture ends up aromatic
NONE = 'none'    # no bond moves


def expect(atoms, bonds, want):
    """Run `thiele` and return `(result, got, wanted)` in fixture indices."""
    mol, ids = build(atoms, bonds)
    result = thiele(mol)
    if want is ALL:
        wanted = {frozenset((i, j)) for i, j, _ in bonds}
    elif want is NONE:
        wanted = set()
    else:
        wanted = set().union(*(cyc(seq) for seq in want))
    return mol, ids, result, aromatic(mol, ids, bonds), wanted


# ------------------------------------------------------------------------------------------------
# THE MOLECULES.  Every one of them is a public compound.
BENZENE = ([('C', 1)] * 6, ring(6, '=-=-=-'))
# N last so that the heteroatom is not also atom 0 in every fixture
PYRIDINE = ([('C', 1)] * 5 + [('N', 0)], ring(6, '=-=-=-'))
PYRIDINIUM = ([('C', 1)] * 5 + [('N', 1, 1)], ring(6, '=-=-=-'))
# a five-ring donor at index 0 and four sp2 carbons: one skeleton, six elements
PYRROLE = ([('N', 1)] + [('C', 1)] * 4, ring(5, '-=-=-'))
FURAN = ([('O', 0)] + [('C', 1)] * 4, ring(5, '-=-=-'))
THIOPHENE = ([('S', 0)] + [('C', 1)] * 4, ring(5, '-=-=-'))
SELENOPHENE = ([('Se', 0)] + [('C', 1)] * 4, ring(5, '-=-=-'))
TELLUROPHENE = ([('Te', 0)] + [('C', 1)] * 4, ring(5, '-=-=-'))
PHOSPHOLE = ([('P', 1)] + [('C', 1)] * 4, ring(5, '-=-=-'))
CYCLOPENTADIENIDE = ([('C', 1, -1)] + [('C', 1)] * 4, ring(5, '-=-=-'))
# the same skeleton with the sign flipped: 4 pi instead of 6, and the file's cleanest refusal
CYCLOPENTADIENYL_CATION = ([('C', 1, 1)] + [('C', 1)] * 4, ring(5, '-=-=-'))
# neutral boron is the empty-orbital atom, so borole is 4 pi where borepine is 6
BOROLE = ([('B', 1)] + [('C', 1)] * 4, ring(5, '-=-=-'))
BOREPINE = ([('B', 1)] + [('C', 1)] * 6, ring(7, '-=-=-=-'))
TROPYLIUM = ([('C', 1, 1)] + [('C', 1)] * 6, ring(7, '-=-=-=-'))
# an NH donor in a seven-ring: 8 pi, the negative that proves Huckel is really being applied
AZEPINE = ([('N', 1)] + [('C', 1)] * 6, ring(7, '-=-=-=-'))
# 1H-imidazole: NH at 0, the other N at 2 carrying a ring double
IMIDAZOLE = ([('N', 1), ('C', 1), ('N', 0), ('C', 1), ('C', 1)], ring(5, '-=-=-'))
# 1H-pyrazole: the two nitrogens adjacent, one donor and one MUST
PYRAZOLE = ([('N', 1), ('N', 0), ('C', 1), ('C', 1), ('C', 1)], ring(5, '-=-=-'))
OXAZOLE = ([('O', 0), ('C', 1), ('N', 0), ('C', 1), ('C', 1)], ring(5, '-=-=-'))
PYRIMIDINE = ([('C', 1), ('N', 0), ('C', 1), ('N', 0), ('C', 1), ('C', 1)], ring(6, '=-=-=-'))
# pyridine N-oxide: three neighbours on the nitrogen and an exocyclic bond that is SINGLE, which is
# what lets it through the exocyclic-double filter that stops every quinone
PYRIDINE_N_OXIDE = ([('C', 1)] * 5 + [('N', 0, 1), ('O', 0, -1)],
                    ring(6, '=-=-=-') + [(5, 6, 1)])

# --- fused and joined.  0 and 1 are the shared atoms wherever two rings meet.
NAPHTHALENE = ([('C', 0)] * 2 + [('C', 1)] * 8,
               [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                (1, 6, 1), (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 0, 1)])
# a five-ring fused to a seven-ring, neither of them 4n+2 on its own: the molecule that proves
# Huckel cannot be applied per ring
AZULENE = ([('C', 0)] * 2 + [('C', 1)] * 8,
           [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 0, 1),
            (1, 5, 1), (5, 6, 2), (6, 7, 1), (7, 8, 2), (8, 9, 1), (9, 0, 2)])
INDOLE = ([('C', 0)] * 2 + [('C', 1)] * 4 + [('N', 1), ('C', 1), ('C', 1)],
          [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1), (5, 0, 2),
           (1, 6, 1), (6, 7, 1), (7, 8, 2), (8, 0, 1)])
# pyrrolo[1,2-a]imidazole, `N1C=CN2C=CC=C12`: TWO donor nitrogens in one five-ring, which a per-ring
# pre-filter caps at one.  Eight atoms share ten pi electrons -- three ring double bonds and a lone
# pair from each nitrogen -- so the whole bicycle is aromatic, as indolizine's is.  Per-ring
# bookkeeping cannot see it: the bridging nitrogen (index 3) spends its pair on the system rather than
# on either ring, and the fusion carbon (index 7) pairs with a carbon in the OTHER ring, so the
# imidazole ring reads as two donors and one double bond on its own.
PYRROLOIMIDAZOLE = ([('N', 1), ('C', 1), ('C', 1), ('N', 0), ('C', 1), ('C', 1), ('C', 1), ('C', 0)],
                    [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 1), (6, 7, 2),
                     (7, 0, 1), (7, 3, 1)])
# tetralin: a benzene fused to a ring with four sp3 carbons.  The saturated ring must not poison the
# aromatic one, which is what the PER-RING pre-filter is for
TETRALIN = ([('C', 0)] * 2 + [('C', 1)] * 4 + [('C', 2)] * 4,
            [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
             (1, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 0, 1)])
# 1,4-naphthoquinone: the same shape with the second ring carrying two exocyclic C=O
NAPHTHOQUINONE = ([('C', 0)] * 2 + [('C', 1)] * 4 + [('C', 0), ('C', 1), ('C', 1), ('C', 0),
                                                     ('O', 0), ('O', 0)],
                  [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                   (1, 6, 1), (6, 7, 1), (7, 8, 2), (8, 9, 1), (9, 0, 1),
                   (6, 10, 2), (9, 11, 2)])
# biphenyl: two systems, one molecule, and the bond between them must stay single
BIPHENYL = ([('C', 0)] + [('C', 1)] * 5 + [('C', 0)] + [('C', 1)] * 5,
            ring(6, '=-=-=-') + ring(6, '=-=-=-', offset=6) + [(0, 6, 1)])

# --- mixtures.  Two components in one container, which is the only way a REFUSED system can be a
# proper subset of the molecule -- and that is what gives the refusal sweep below any resolution.
CATION_AND_ETHENE = (CYCLOPENTADIENYL_CATION[0] + [('C', 2), ('C', 2)],
                     CYCLOPENTADIENYL_CATION[1] + [(5, 6, 2)])
BENZENE_AND_CATION = (BENZENE[0] + CYCLOPENTADIENYL_CATION[0],
                      BENZENE[1] + [(6 + i, 6 + j, o) for i, j, o in CYCLOPENTADIENYL_CATION[1]])

# --- never candidates.
CYCLOBUTADIENE = ([('C', 1)] * 4, ring(4, '=-=-'))
CYCLOOCTATETRAENE = ([('C', 1)] * 8, ring(8, '=-=-=-=-'))
CYCLOHEXENE = ([('C', 1)] * 2 + [('C', 2)] * 4, ring(6, '=-----'))
CYCLOHEXANE = ([('C', 2)] * 6, ring(6, '------'))
# fulvene: every ring atom carries exactly one double bond and one of them points OUT
FULVENE = ([('C', 0)] + [('C', 1)] * 4 + [('C', 2)], ring(5, '-=-=-') + [(0, 5, 2)])
# 2-pyridone and 4-pyridone: an NH donor, a perfect Kekule match, and an exocyclic carbonyl
PYRIDONE_2 = ([('N', 1), ('C', 0), ('C', 1), ('C', 1), ('C', 1), ('C', 1), ('O', 0)],
              ring(6, '--=-=-') + [(1, 6, 2)])
PYRIDONE_4 = ([('C', 0), ('C', 1), ('C', 1), ('N', 1), ('C', 1), ('C', 1), ('O', 0)],
              ring(6, '-=--=-') + [(0, 6, 2)])
P_BENZOQUINONE = ([('C', 0), ('C', 1), ('C', 1), ('C', 0), ('C', 1), ('C', 1),
                   ('O', 0), ('O', 0)],
                  ring(6, '-=--=-') + [(0, 6, 2), (3, 7, 2)])
# borazine, all single bonds: no Kekule form to match at all.  IT IS THE ONE NEGATIVE THE MATCHING
# CHECK CANNOT SEE -- three N donors want no ring double bond and find none, three neutral borons
# contribute an empty orbital, and six pi over six atoms passes Huckel.  What refuses it is the
# pre-filter's requirement that a candidate ring hold at least one double bond of its own, which is
# the same sentence as "the caller is holding a Kekule form".
BORAZINE = ([('B', 1), ('N', 1)] * 3, ring(6, '------'))
# 1,4-dihydropyridine and 2H-pyran: a donor heteroatom AND an sp3 carbon.  What refuses them is not a
# count of non-sp2 atoms -- PYRROLOIMIDAZOLE above has two and is aromatic -- but the sp3 carbon
# itself: a neutral carbon with two hydrogens has no valid aromatic state at all.
DIHYDROPYRIDINE = ([('N', 1), ('C', 1), ('C', 1), ('C', 2), ('C', 1), ('C', 1)],
                   ring(6, '-=--=-'))
PYRAN_2H = ([('O', 0), ('C', 2), ('C', 1), ('C', 1), ('C', 1), ('C', 1)], ring(6, '--=-=-'))


POSITIVE = [('benzene', BENZENE, ALL),
            ('pyridine', PYRIDINE, ALL),
            ('pyridinium', PYRIDINIUM, ALL),
            ('pyrrole', PYRROLE, ALL),
            ('furan', FURAN, ALL),
            ('thiophene', THIOPHENE, ALL),
            ('selenophene', SELENOPHENE, ALL),
            ('tellurophene', TELLUROPHENE, ALL),
            ('phosphole', PHOSPHOLE, ALL),
            ('cyclopentadienide', CYCLOPENTADIENIDE, ALL),
            ('borepine', BOREPINE, ALL),
            ('tropylium', TROPYLIUM, ALL),
            ('imidazole', IMIDAZOLE, ALL),
            ('pyrazole', PYRAZOLE, ALL),
            ('oxazole', OXAZOLE, ALL),
            ('pyrimidine', PYRIMIDINE, ALL),
            ('pyridine N-oxide', PYRIDINE_N_OXIDE, [range(6)]),
            ('naphthalene', NAPHTHALENE, ALL),
            ('azulene', AZULENE, ALL),
            ('indole', INDOLE, ALL),
            ('pyrrolo[1,2-a]imidazole', PYRROLOIMIDAZOLE, ALL),
            ('tetralin', TETRALIN, [range(6)]),
            ('naphthoquinone', NAPHTHOQUINONE, [range(6)]),
            ('biphenyl', BIPHENYL, [range(6), range(6, 12)])]

NEGATIVE = [('cyclobutadiene', CYCLOBUTADIENE),
            ('cyclooctatetraene', CYCLOOCTATETRAENE),
            ('cyclohexene', CYCLOHEXENE),
            ('cyclohexane', CYCLOHEXANE),
            ('fulvene', FULVENE),
            ('2-pyridone', PYRIDONE_2),
            ('4-pyridone', PYRIDONE_4),
            ('p-benzoquinone', P_BENZOQUINONE),
            ('borazine', BORAZINE),
            ('1,4-dihydropyridine', DIHYDROPYRIDINE),
            ('2H-pyran', PYRAN_2H),
            ('cyclopentadienyl cation', CYCLOPENTADIENYL_CATION),
            ('borole', BOROLE),
            ('1H-azepine', AZEPINE)]

# the three negatives that WERE proposed and declined, with the ring that was declined
DECLINED = [('cyclopentadienyl cation', CYCLOPENTADIENYL_CATION, frozenset(range(5))),
            ('borole', BOROLE, frozenset(range(5))),
            ('1H-azepine', AZEPINE, frozenset(range(7)))]


# ------------------------------------------------------------------------------------------------
# WHAT BECOMES AROMATIC.
@mark.parametrize('name,fixture,want', POSITIVE)
def test_the_aromatic_bonds_are_exactly_these(name, fixture, want):
    """One row per molecule, and the assertion is on the WHOLE bond set rather than on a count.

    A count would pass for a mechanism that aromatised the wrong six bonds, and three of these rows
    exist precisely because a bond next door must NOT move: pyridine N-oxide's N-O, tetralin's and
    naphthoquinone's saturated ring, biphenyl's inter-ring bond.
    """
    mol, ids, result, got, wanted = expect(*fixture, want)
    assert got == wanted
    assert result.changed
    assert notices(result) == [] and result.refused == []


@mark.parametrize('name,fixture', NEGATIVE)
def test_nothing_moves(name, fixture):
    """The other half.  `changed` is False and no bond of the fixture holds order 4."""
    mol, ids, result, got, wanted = expect(*fixture, NONE)
    assert got == set()
    assert not result.changed


@mark.parametrize('name,fixture,expected', DECLINED)
def test_a_declined_system_is_named_and_logged(name, fixture, expected):
    """A Kekule match, every atom satisfied, and the wrong electron count -- so it is REPORTED.

    This is the distinction the file's docstring opens with.  Silence would be indistinguishable
    from "your ring was never a candidate", and these three are the cases where the caller drew
    something that looks exactly like an aromatic ring.
    """
    mol, ids, result, got, _ = expect(*fixture, NONE)
    assert got == set()
    assert refused(mol, ids, result) == {expected}
    assert len(notices(result)) == 1
    assert 'not a Kekule form' in notices(result)[0]


@mark.parametrize('name,fixture', [row for row in NEGATIVE
                                   if row[0] not in {n for n, _, _ in DECLINED}])
def test_a_ring_that_was_never_a_candidate_is_not_reported(name, fixture):
    """The complement of the test above, and the reason `.refused` is worth reading.

    Every ring here is dropped before any electron is counted -- by size, by an sp3 atom, by an
    exocyclic double bond, or by having no double bonds at all.  Reporting them would put a line in
    `.log` for every cyclohexane in a corpus.
    """
    mol, ids, result, _, _ = expect(*fixture, NONE)
    assert result.refused == [] and notices(result) == []


def test_the_PRE_FILTER_names_the_two_states_it_declines_a_good_ring_for():
    """A benzene with one unknown hydrogen count must not come back `(False, [], [])`.

    The pre-filter runs before any component exists, so a ring it drops is never reached by the loop
    that fills `.refused` -- and two of its five exits are refusals in the docstring's sense: the ring
    IS a Kekule aromatic and is declined for a state that a caller can fix.  Returning either silently
    breaks the acceptance bar in the one place this project says a refusal must name itself.

    THE UNKNOWN COUNT IS THE CASE THAT MATTERS, and it is ordinary rather than exotic: between a read
    and a `kekule()` the pyrrole-versus-pyridine atom is exactly this, and a molecule built atom by
    atom in an edit session has every count unknown until `derive_hydrogens()` runs.  Without the
    diagnostic, such a molecule aromatises to nothing and says nothing about the counts being the
    reason rather than the bond orders.

    Benzene is the subject for both halves precisely because it is unimpeachable: whatever else the
    filter thinks, six carbons with alternating orders are a Kekule aromatic ring, so a `False` here
    can only be the state under test.
    """
    atoms = [('C', 1)] * 6
    bonds = ring(6, '=-=-=-')

    mol, ids = build(atoms, bonds)
    with mol.edit() as e:
        e.set_hydrogens(ids[0], H_UNKNOWN)
    result = thiele(mol)
    assert not result.changed
    assert refused(mol, ids, result) == {frozenset(range(6))}, 'the ring, by stable id'
    assert len(notices(result)) == 1
    assert 'unknown implicit hydrogen count' in notices(result)[0]
    assert f'atom {ids[0]}' in notices(result)[0], 'the offending atom is named, not just the ring'
    assert 'derive_hydrogens()' in notices(result)[0], 'and what to do about it'

    mol, ids = build(atoms, bonds)
    with mol.edit() as e:
        e.set_radical(ids[0], True)
    result = thiele(mol)
    assert not result.changed
    assert refused(mol, ids, result) == {frozenset(range(6))}
    assert len(notices(result)) == 1
    assert 'radical' in notices(result)[0] and f'atom {ids[0]}' in notices(result)[0]

    # ...and the same benzene with neither state is still silent, so the two assertions above are
    # about the state and not about benzene having acquired a log line
    mol, ids = build(atoms, bonds)
    result = thiele(mol)
    assert result.changed and result.refused == [] and notices(result) == []


# ------------------------------------------------------------------------------------------------
# THE TWO MOLECULES THAT PIN THE TWO DESIGN DECISIONS.  Both are already rows above; both get a
# named test as well, because a row that silently starts passing for the wrong reason is invisible.
def test_tetralin_keeps_its_benzene_although_the_fused_ring_is_saturated():
    """The PER-RING pre-filter, which is what makes this work at all.

    Without it the saturated ring's atoms would join the candidate support, and then the four sp3
    carbons would fail the Kekule match for the whole component -- one bad fused ring would cost a
    perfectly good benzene its aromaticity.
    """
    mol, ids, result, got, _ = expect(*TETRALIN, [range(6)])
    assert got == cyc(range(6))
    assert mol.order_of(ids[6], ids[7]) == 1
    assert mol.order_of(ids[1], ids[6]) == 1     # the bond joining the two rings


def test_azulene_is_aromatic_although_neither_of_its_rings_is_4n_plus_2():
    """5 pi and 7 pi separately, 10 pi together -- and 10 % 4 == 2.

    Huckel is applied to a candidate system that is a SINGLE CYCLE, and azulene's is not: the two
    shared atoms have candidate degree 3, so the rule is skipped and the Kekule match decides.  A
    per-ring Huckel test would refuse this molecule, and it is aromatic.
    """
    mol, ids, result, got, wanted = expect(*AZULENE, ALL)
    assert got == wanted and result.changed


def test_every_system_that_comes_out_says_so():
    """The aromatisation itself is the event `mol.log` exists to hold, so it is a record.

    ONE PER ACCEPTED SYSTEM, `INFO`, naming that system's atoms in sorted stable ids and the electron
    count that admitted it -- biphenyl is two lines and naphthalene is one, which is the same partition
    the aromatic bond set is asserted against above.  A refused system gets its refusal and nothing
    else, and a molecule already aromatic says nothing at all.
    """
    mol, ids = build(*BIPHENYL)
    result = thiele(mol)
    assert [(x.rule, x.atoms, x.severity) for x in result.log] == \
        [('thiele:aromatized', tuple(ids[i] for i in range(6)), 'info'),
         ('thiele:aromatized', tuple(ids[i] for i in range(6, 12)), 'info')]
    assert all('6 pi electrons' in x for x in result.log)
    second = thiele(mol)
    assert not second.changed and second.log == [], 'an aromatic molecule said something anyway'

    mol, ids = build(*NAPHTHALENE)
    result = thiele(mol)
    assert [x.atoms for x in result.log] == [tuple(ids[i] for i in range(10))]
    assert '10 pi electrons' in result.log[0], 'the count Huckel was applied to is part of the record'

    # a refused system beside an accepted one: one line each, and only the benzene's is the rewrite
    mol, ids = build(*BENZENE_AND_CATION)
    result = thiele(mol)
    assert [x.rule for x in result.log] == ['thiele:aromatized', 'thiele:not-kekule-form']
    assert [x.atoms for x in result.log if x.rule == 'thiele:aromatized'] == \
        [tuple(ids[i] for i in range(6))]


def test_borole_is_refused_where_borepine_is_accepted():
    """Neutral boron brings an EMPTY orbital, so a boron ring counts as one atom short of a donor.

    Five-ring: 4 pi, refused.  Seven-ring: 6 pi, accepted.  chython 2 aromatises borole -- its filter
    does not count electrons -- so a template ported from V2 changes answer on this ring.
    """
    _, _, refusal, borole, _ = expect(*BOROLE, NONE)
    _, _, accept, borepine, wanted = expect(*BOREPINE, ALL)
    assert borole == set() and not refusal.changed
    assert borepine == wanted and accept.changed


# ------------------------------------------------------------------------------------------------
# THE METHOD.  One operation, one name, two spellings -- and the test is that they are one operation.
def test_the_method_is_the_function_and_not_a_second_implementation():
    """`mol.thiele()` and `thiele(mol)` are the same call, which is why only one of them is tested.

    Every other test in this file goes through the function; the method is what callers actually
    write.  Asserting the equivalence once here is what makes those tests cover both -- and the
    failure it guards is the ordinary one, a method that grows a defaulted argument or a pre-check the
    function does not have, at which point every test above stops describing the thing being used.
    """
    mol, ids = build(*BENZENE)
    assert aromatic(mol, ids, BENZENE[1]) == set(), 'the fixture starts Kekule'
    result = mol.thiele()
    assert isinstance(result, ThieleResult)
    assert result.changed
    assert notices(result) == [] and result.refused == []
    assert aromatic(mol, ids, BENZENE[1]) == cyc(range(6))
    # and the pair is symmetric at the method surface too, since a caller who found one of them by
    # autocomplete has to find the other the same way
    assert mol.kekule().changed
    assert aromatic(mol, ids, BENZENE[1]) == set()


# ------------------------------------------------------------------------------------------------
# IDEMPOTENCE AND THE ROUND TRIP.
def stated(atoms, ids):
    """The fixture's hydrogen counts, keyed by stable id, as `kekule` wants them."""
    return {ids[k]: atoms[k][1] for k in range(len(atoms))}


@mark.parametrize('name,fixture,want', POSITIVE)
def test_running_it_twice_moves_nothing_the_second_time(name, fixture, want):
    """`changed` is measured against the arena, so the second call reports False rather than
    promising it.

    The mechanism is worth naming: an order-4 bond fails the per-ring pre-filter, so an
    already-aromatic ring is not a candidate at all.  Idempotence here is not a special case in the
    code, it is the same rule refusing to aromatise what is already aromatic.
    """
    mol, ids, first, got, wanted = expect(*fixture, want)
    again = thiele(mol)
    assert not again.changed
    assert again.log == [] and again.refused == []
    assert aromatic(mol, ids, fixture[1]) == got


@mark.parametrize('name,fixture,want', POSITIVE)
def test_the_round_trip_through_kekule_is_the_identity_on_the_aromatic_set(name, fixture, want):
    """`thiele` then `kekule` then `thiele` finds the same bonds aromatic.

    Not the same BOND ORDERS in between -- which Kekule form comes back is a free choice between
    equivalent answers -- so the invariant is stated on the aromatic edge set, which is the thing
    both operations agree about.  The hydrogen counts are handed to `kekule` explicitly; the test
    below is why.
    """
    mol, ids, _, first, wanted = expect(*fixture, want)
    assert first == wanted
    bonds = fixture[1]
    edges = [(ids[i], ids[j]) for i, j, _ in bonds if frozenset((i, j)) in first]
    assert notices(kekule(mol, edges, stated(fixture[0], ids))) == []
    assert aromatic(mol, ids, bonds) == set()
    assert thiele(mol).changed
    assert aromatic(mol, ids, bonds) == first


@mark.parametrize('name,fixture', [('imidazole', IMIDAZOLE), ('pyrazole', PYRAZOLE)])
def test_kekules_no_argument_form_ROUND_TRIPS_a_two_nitrogen_five_ring(name, fixture):
    """`kekule()` with no `stated_h` reads the STORED counts, which is what lets these two survive it.

    Two things carry that, and only the second is `_kekule.pxi`'s.  `add_atom(implicit_h=None)` stores
    `H_UNKNOWN` rather than 0, so "nobody said" is a value a reader can recognise; and `arom_setup`
    takes the stored nibble as its DEFAULT SOURCE, with `stated_h` overriding it, rather than forcing
    every atom absent from the dict to unstated.  A builder mid-flight is still treated as unstated,
    since every nibble it has not written holds `H_UNKNOWN`.

    Read `stated_h=None` as "nobody stated any count" instead and every two-coordinate aromatic N
    becomes the free choice; on a ring with TWO of them the freedom picks the wrong one -- imidazole
    comes back with a double bond on the N holding the hydrogen (a four-valent neutral N), pyrazole
    with neither N doubled.  `thiele` then declines the ring and the molecule is silently no longer
    imidazole.

    The H counts are asserted UNCHANGED across the round trip because `kekule` reads them: a version
    that WROTE them would also pass the ring assertions.
    """
    atoms, bonds = fixture
    mol, ids = build(atoms, bonds)
    assert thiele(mol).changed
    before = aromatic(mol, ids, bonds)
    assert before != set()
    assert kekule(mol).changed
    assert [mol.implicit_h_of(ids[k]) for k in range(len(atoms))] == [a[1] for a in atoms], \
        'kekule reads the stored counts and must not rewrite them'
    again = thiele(mol)
    assert again.changed, 'the kekule form is a real one, so thiele takes the ring back'
    assert aromatic(mol, ids, bonds) == before, 'and takes back the SAME ring, bond for bond'


# ------------------------------------------------------------------------------------------------
# THE EDGES OF THE ENTRY POINT.
def test_an_acyclic_molecule_is_a_no_op():
    """No ring segment to read, and the early return says so without allocating scratch."""
    mol, ids = build([('C', 3), ('C', 2), ('O', 1)], [(0, 1, 1), (1, 2, 1)])
    result = thiele(mol)
    assert not result.changed and notices(result) == [] and result.refused == []


def test_an_empty_molecule_is_a_no_op():
    result = thiele(MoleculeContainer())
    assert not result.changed and notices(result) == [] and result.refused == []


def test_pending_edits_are_refused():
    """Same contract as `kekule`: the operation reads the built structure, so the journal must be
    clean before it runs."""
    mol, ids = build(*BENZENE)
    with raises(RuntimeError):
        with mol.edit():
            mol.add_atom('C')
            thiele(mol)


def test_the_result_is_a_named_object_and_not_a_tuple():
    """A caller unpacking three fields positionally is the call site that breaks when a fourth
    arrives, so unpacking must not work at all."""
    mol, ids = build(*BENZENE)
    result = thiele(mol)
    assert isinstance(result, ThieleResult)
    assert not isinstance(result, tuple)
    with raises(TypeError):
        a, b, c = result
    assert 'changed=True' in repr(result)


# ------------------------------------------------------------------------------------------------
# INVARIANCE.  The verdict is a property of the molecule, so it must not depend on the order the
# atoms were created in -- which decides slot numbers, CSR layout, ring-perception order and the
# order the components are visited in.
def creation_orders(n, cap=5040, samples=200, seed=0):
    """Every permutation while that is affordable, else a fixed pseudo-random sample of them."""
    if factorial(n) <= cap:
        return list(permutations(range(n)))
    rnd = Random(seed)
    out = []
    for _ in range(samples):
        order = list(range(n))
        rnd.shuffle(order)
        out.append(tuple(order))
    return out


def sweep(fixture, key='index'):
    """The set of answers over the creation orders: `(aromatic edges, refused systems)`.

    `key='index'` reads both back in fixture indices, which is the molecule's own frame and must give
    ONE answer.  `key='sid'` reads them in stable ids, which is the creation order's frame -- that is
    the same trick `test_smiles_write_h_unknown.py` plays with the `i` spec, and it exists to prove
    the sweep is capable of reporting more than one thing.
    """
    atoms, bonds = fixture
    answers = set()
    for order in creation_orders(len(atoms)):
        mol, ids = build(atoms, bonds, order)
        result = thiele(mol)
        if key == 'index':
            answers.add((frozenset(aromatic(mol, ids, bonds)),
                         frozenset(refused(mol, ids, result))))
        else:
            answers.add((frozenset(frozenset((ids[i], ids[j])) for i, j, _ in bonds
                                   if mol.order_of(ids[i], ids[j]) == 4),
                         frozenset(frozenset(system) for system in result.refused)))
    return answers


# (name, fixture) -- one monocycle, one heterocycle, one refusal, and the fused cases where the
# component walk and the ring pre-filter interact
SWEEPS = [('benzene', BENZENE),
          ('pyrrole', PYRROLE),
          ('imidazole', IMIDAZOLE),
          ('cyclopentadienyl cation', CYCLOPENTADIENYL_CATION),
          ('cation and ethene', CATION_AND_ETHENE),
          ('benzene and cation', BENZENE_AND_CATION),
          ('azulene', AZULENE),
          ('tetralin', TETRALIN),
          ('naphthoquinone', NAPHTHOQUINONE),
          ('biphenyl', BIPHENYL)]


@mark.parametrize('name,fixture', SWEEPS)
def test_the_verdict_is_the_same_from_every_creation_order(name, fixture):
    """One answer, and `refused` is swept along with the aromatic set.

    The refusal row matters as much as the acceptances: a Huckel count that came out order-dependent
    would mean the component walk was visiting a different set of atoms, and the accepted molecules
    could not show that -- they would simply all be aromatic.
    """
    assert len(sweep(fixture)) == 1


# the MEASURED number of distinct answers in the creation order's own frame, per fixture.  Every
# number here was measured and none was predicted -- the first guesses were 24 and 60 for the
# five-rings and both were wrong, because what varies is the labelling of a CYCLE and there are
# 5! / (2 * 5) = 12 of those, not 5! / 5
CAN_FAIL = [('benzene', BENZENE, 60),
            ('pyrrole', PYRROLE, 12),
            ('imidazole', IMIDAZOLE, 12),
            ('cation and ethene', CATION_AND_ETHENE, 21)]


@mark.parametrize('name,fixture,n', CAN_FAIL)
def test_the_sweep_can_fail(name, fixture, n):
    """Ruling F102: the sweep above is evidence only if it could have reported more than one answer.

    Read in stable ids the same molecule gives many answers -- 60 for benzene, which is 720 creation
    orders divided by the 12 automorphisms of a six-cycle, 12 for each five-ring, and 21 for the
    mixture, which is the number of five-atom subsets of seven that its refused ring can land on.
    Those numbers are the sweep's resolution: a mechanism that answered the same thing in every frame
    would give 1 here too, and then the test above would be measuring nothing.

    The mixture is in this list and the bare cation is not, for the reason the next test records.
    """
    assert len(sweep(fixture, 'sid')) == n


def test_the_refusal_sweep_has_no_resolution_on_a_single_component_molecule():
    """MEASURED, and the reason `CATION_AND_ETHENE` exists.

    A refused system is a set of atoms, and when it is ALL of them the set is the same in every frame
    -- so sweeping the bare cation in stable ids gives one answer and proves nothing about the
    refusal.  Adding an ethene the mechanism does not care about makes the refused ring a proper
    subset, and the sweep can then tell frames apart.  Recording the degenerate measurement rather
    than deleting the row, because "the sweep passed" reads the same either way and this is the
    difference between a test and a decoration.
    """
    assert len(sweep(CYCLOPENTADIENYL_CATION, 'sid')) == 1
    assert len(sweep(CATION_AND_ETHENE, 'sid')) > 1


def test_a_mixture_decides_component_by_component():
    """Benzene and the cyclopentadienyl cation in one container: one aromatised, one refused.

    The per-component loop is what this asserts -- not "some bonds moved", but that the accepted
    component moved ALL of its bonds and the refused one none of its own, in the same call, with the
    refusal reported.  A mechanism that gave up on the whole molecule after one refusal would pass
    every other test in this file.
    """
    atoms, bonds = BENZENE_AND_CATION
    mol, ids = build(atoms, bonds)
    result = thiele(mol)
    assert result.changed
    assert aromatic(mol, ids, bonds) == cyc(range(6))
    assert refused(mol, ids, result) == {frozenset(range(6, 11))}
    assert len(notices(result)) == 1
