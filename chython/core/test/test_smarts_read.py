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
"""The SMARTS reader: every primitive, every refusal, and the places chython 2 read it differently.

Three kinds of test live here and they are not interchangeable.

**Primitives.**  One spelling, one panel of public molecules, one expected hit set.  These pin the
*meaning* of a letter, so they are written against molecules and not against the journal: a test that
asserted "this emitted PRIM_DEGREE" would pass just as happily if PRIM_DEGREE meant something else.

**Refusals.**  A malformed string raises `IncorrectSmarts` and the message ends in a byte offset.  The
tests assert the offset, because "somewhere in your 400-character template" is not a diagnostic.
A string that is merely *odd* is not in this section -- it is parsed, and if anything was dropped a
line says so in the log.  That split is the reader's whole error policy: syntax raises, chemistry is
logged, and nothing is silently repaired.

**Acceptance tests for the V2 divergences.**  A V2 reading is a specification input, and the test that
pins the V3 reading is where it goes.  Each one below names what V2 does, MEASURED against chython 2.24
and not recalled:

* `~` matches nothing in V2 (it lexes to bond order 8, the dative bond).  Here it is a real five-way
  disjunction.
* `c1ccccc1` in V2 matches CYCLOHEXANE and not benzene: its tokeniser collapses the aromatic-atom
  token onto the plain one, dropping the aromatic demand and leaving the single bonds behind.
* `[C;H3]`, `[c]`, `[*]`, `*`, `[C;R]`, `[C;r0]`, `[!C]`, `[C&D1]`, `[CD1]` and `z5`/`z6` are all
  parse errors in V2, which has five primitives and one acyclic spelling (`D h r !R a`), no negated
  element, no `&` and no juxtaposition.
* V2 reads a contradiction as its last clause -- `[C;D1;D2]` there IS `[C;D2]`, measured -- and
  accepts a map number past what any arena can store.  Both are refused here, at the position where
  they were written.
* `;M` in V2 hands the atom a map number from a module-level counter.  Masking is a flag here.

And one divergence that is neither side's error: `z`.  V2 saturates -- sp3, promoted once per
double bond, capped at 3 -- so a sulfone S, a nitro N and an allene's middle carbon all answered `z3`
there beside genuine sp carbons.  The core reports what it found, `5` for two cumulated doubles and
`6` for anything past that, and `z3` means sp and nothing else.  The reader passes the number through
unchanged; translating the 50 `z3` primitives carried over from V2 is a separate task with its own
written mapping, and this file only makes the difference visible.

The last section runs the whole thing against chython 2.24 in a separate interpreter -- the ORACLE,
pinned to that version and skipped cleanly when it is not installed.  That subprocess runs under `-I`
and asserts on the module path it actually imported, because a child started with plain `-c` from the
repo root imports the tree under test: the differential then compares this reader against itself and
passes with flying colours, which is the one failure mode a differential cannot survive.  Two
exclusions apply there and
both are the matching kernel's semantics rather than the reader's, verified by building the same query
through the journal API by hand where no lexer is involved:

1. The core's `is_substructure` is a MONOMORPHISM, as Daylight SMARTS specifies -- `CCC` matches
   cyclopropane.  V2's is induced: it rejects a match whose mapped atoms carry a bond the query did
   not state.
2. `.` here means "no bond stated", again as Daylight specifies.  V2 additionally demands that the
   two sides land in different connected components, so its `[O;D1].[O;D1]` misses benzoic acid.
"""
from json import dumps, loads
from os import environ
from pathlib import Path
from re import compile as compile_regex

from pytest import mark, raises, skip

from chython.core import IncorrectSmarts, QueryContainer, read_smarts, read_smiles
from chython.core.test.oracle import ask


# The panel.  Public compounds and minimal probes; every name below is what the molecule IS, so a
# failing parametrised case reads as a sentence.
PANEL = {
    'methane': 'C', 'water': 'O', 'iron': '[Fe]',
    'propane': 'CCC', 'propene': 'CC=C', 'propyne': 'CC#C', 'allene': 'C=C=C',
    'ethanol': 'CCO', 'acetone': 'CC(=O)C', 'acetic_acid': 'CC(=O)O',
    'acetonitrile': 'CC#N', 'ethylamine': 'CCN', 'bromoethane': 'CCBr',
    'dimethyl_sulfone': 'CS(=O)(=O)C', 'dimethyl_sulfide': 'CSC',
    'benzene': 'c1ccccc1', 'toluene': 'Cc1ccccc1', 'pyridine': 'c1ccncc1',
    'pyrrole': 'c1cc[nH]c1', 'phenol': 'Oc1ccccc1', 'aniline': 'Nc1ccccc1',
    'chlorobenzene': 'Clc1ccccc1', 'naphthalene': 'c1ccc2ccccc2c1',
    'nitrobenzene': 'O=[N+]([O-])c1ccccc1', 'benzoic_acid': 'OC(=O)c1ccccc1',
    'cyclopropane': 'C1CC1', 'cyclopentane': 'C1CCCC1', 'cyclohexane': 'C1CCCCC1',
    'thf': 'C1CCOC1', 'deuteriomethane': '[2H]C', 'carbon13_methane': '[13CH4]',
}
MOLECULES = {k: read_smiles(v) for k, v in PANEL.items()}

# The second panel, for the shipped-template differential at the bottom of the file only.  The panel
# above is built of minimal probes, which is what a test about one primitive needs and exactly the
# wrong thing for a test about real templates: against it the 486 templates this tree ships match 59
# cells out of 15066, so the differential mostly agrees that neither reader matched -- which two
# readers can do while disagreeing about everything.  These are
# textbook reagents and bench chemicals chosen so the shipped functional-group patterns actually
# fire.  Public compounds only.
REAGENTS = {
    'boc_glycine': 'OC(=O)CNC(=O)OC(C)(C)C', 'boc_piperidine': 'O=C(OC(C)(C)C)N1CCCCC1',
    'boc_anhydride': 'O=C(OC(C)(C)C)OC(=O)OC(C)(C)C',
    'cbz_alanine': 'CC(NC(=O)OCc1ccccc1)C(=O)O', 'benzyl_alcohol': 'OCc1ccccc1',
    'phenylboronic_acid': 'OB(O)c1ccccc1', 'phenylboronic_pinacol_ester': 'CC1(C)OB(c2ccccc2)OC1(C)C',
    'bromobenzene': 'Brc1ccccc1', 'iodobenzene': 'Ic1ccccc1', 'benzyl_bromide': 'BrCc1ccccc1',
    'bromopyridine': 'Brc1ccncc1', 'chloropyridine': 'Clc1ccccn1', 'bromobutane': 'CCCCBr',
    'chloronitrobenzene': 'Clc1ccc(cc1)[N+](=O)[O-]',
    'benzamide': 'NC(=O)c1ccccc1', 'acetanilide': 'CC(=O)Nc1ccccc1',
    'ethyl_benzoate': 'CCOC(=O)c1ccccc1', 'ethyl_acetate': 'CCOC(C)=O',
    'ethyl_bromoacetate': 'CCOC(=O)CBr', 'diethyl_malonate': 'CCOC(=O)CC(=O)OCC',
    'ethyl_cyanoacetate': 'CCOC(=O)CC#N',
    'benzaldehyde': 'O=Cc1ccccc1', 'acetophenone': 'CC(=O)c1ccccc1', 'benzonitrile': 'N#Cc1ccccc1',
    'chloroacetone': 'CC(=O)CCl', 'cyclohexanone': 'O=C1CCCCC1', 'cyclohexanol': 'OC1CCCCC1',
    'benzoyl_chloride': 'ClC(=O)c1ccccc1', 'acetic_anhydride': 'CC(=O)OC(C)=O',
    'succinic_anhydride': 'O=C1CCC(=O)O1', 'trifluoroacetic_acid': 'OC(=O)C(F)(F)F',
    'benzenesulfonyl_chloride': 'O=S(=O)(Cl)c1ccccc1', 'benzenesulfonamide': 'NS(=O)(=O)c1ccccc1',
    'tosyl_chloride': 'Cc1ccc(cc1)S(=O)(=O)Cl', 'methyl_tosylate': 'Cc1ccc(cc1)S(=O)(=O)OC',
    'methyl_mesylate': 'COS(C)(=O)=O', 'phenyl_triflate': 'O=S(=O)(Oc1ccccc1)C(F)(F)F',
    'styrene_oxide': 'C1OC1c1ccccc1', 'phenylacetylene': 'C#Cc1ccccc1',
    'phenyl_isocyanate': 'O=C=Nc1ccccc1', 'acetophenone_oxime': 'CC(=NO)c1ccccc1',
    'phenylhydrazine': 'NNc1ccccc1', 'phenylurea': 'NC(=O)Nc1ccccc1', 'thiophenol': 'Sc1ccccc1',
    'anisole': 'COc1ccccc1', 'benzylamine': 'NCc1ccccc1', 'n_methylaniline': 'CNc1ccccc1',
    'glycine': 'NCC(=O)O', 'proline': 'OC(=O)C1CCCN1', 'nicotinic_acid': 'OC(=O)c1cccnc1',
    'imidazole': 'c1c[nH]cn1', 'pyrazole': 'c1cc[nH]n1', 'indole': 'c1ccc2[nH]ccc2c1',
    'thiophene': 'c1ccsc1', 'furan': 'c1ccoc1', 'piperidine': 'C1CCNCC1',
    'morpholine': 'C1COCCN1', 'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'paracetamol': 'CC(=O)Nc1ccc(O)cc1', 'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'caffeine': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
}
TEMPLATE_PANEL = {**PANEL, **REAGENTS}
TEMPLATE_MOLECULES = {k: read_smiles(v) for k, v in TEMPLATE_PANEL.items()}


def hits(pattern):
    """Which panel members the pattern matches, by name and sorted, so a diff reads."""
    q = read_smarts(pattern)
    return sorted(k for k, m in MOLECULES.items() if q.is_substructure(m))


# ----------------------------------------------------------------------------------------------
# ELEMENTS, and the ten letters that are primitives instead
# ----------------------------------------------------------------------------------------------
def test_bare_organic_subset():
    """The organic subset needs no brackets, and a two-letter symbol is read greedily."""
    assert hits('Br') == ['bromoethane']
    assert hits('Cl') == ['chlorobenzene']
    assert read_smarts('CCO').atom_count_sealed() == 3


@mark.parametrize('pattern,expect', [
    ('[C;h2]', 'implicit hydrogens, NOT hydrogen ANDed with something'),
    ('[C;x0]', 'heteroatom count, not xenon'),
    ('[C;a]', 'aromatic, not actinium'),
    ('[C;r6]', 'ring size, not radon'),
    ('[C;z1]', 'hybridization, not zinc'),
    ('[C;D2]', 'degree, not dysprosium'),
    ('[C;R]', 'ring count, not radium'),
    ('[C;H2]', 'total hydrogens, not helium or hydrogen'),
    ('[C;M]', 'masked, not magnesium'),
])
def test_ten_letters_are_primitives_not_elements(pattern, expect):
    """`A D H M R a h r x z` never begin a ONE-letter element symbol inside a bracket.

    Without the rule `[C;h2]` reads as hydrogen and `[C;x0]` as xenon, and both compile -- to a
    query that silently matches nothing, which is the worst of the available outcomes.  Two-letter
    symbols starting with the uppercase five are unaffected; the next test is that half.
    """
    assert read_smarts(pattern).atom_count_sealed() == 1, expect


@mark.parametrize('pattern', ['[Dy]', '[Ho]', '[Mg]', '[Ag]', '[Ru]', '[Ar]', '[At]', '[Hf]'])
def test_two_letter_symbols_still_read(pattern):
    """The primitive-letter rule is about ONE-letter lookups only."""
    assert read_smarts(pattern).atom_count_sealed() == 1


def test_h_a_m_are_the_element_in_first_position():
    """`[H]` is a hydrogen atom, `[A]` any atom, `[M]` any metal -- and only leading the body.

    Later in the same bracket the same letters are `total_h`, nothing, and masked.  The corpus relies
    on the positional reading: `[A:1]` appears in 90 product templates.
    """
    assert hits('[H]') == ['deuteriomethane']            # the only EXPLICIT hydrogen in the panel
    assert hits('[M]') == ['iron']
    assert len(hits('[A]')) == len(PANEL)                # every molecule has some atom
    assert hits('[C;M]') == hits('[C]')                  # masking is not a test
    assert hits('[C;A]') == hits('[C]')                  # `A` after the first primitive says nothing


def test_atomic_number_and_wildcards():
    assert hits('[#7]') == hits('[N]')
    assert len(hits('*')) == len(PANEL)
    assert hits('[*;D4]') == ['dimethyl_sulfone']        # the only four-coordinate heavy atom


def test_element_list_is_one_or_of_elements():
    assert hits('[F,Cl,Br,I]') == ['bromoethane', 'chlorobenzene']
    assert hits('[C,N,O;D1;x0]') == sorted(set(hits('[C;D1;x0]')) | set(hits('[N;D1;x0]'))
                                           | set(hits('[O;D1;x0]')))


def test_negated_element():
    assert 'benzene' not in hits('[!C]')                 # all carbon
    assert 'ethanol' in hits('[!C]')


# ----------------------------------------------------------------------------------------------
# ISOTOPES
# ----------------------------------------------------------------------------------------------
def test_isotope():
    assert hits('[13C]') == ['carbon13_methane']
    assert hits('[2H]') == ['deuteriomethane']
    assert 'carbon13_methane' in hits('[C]')             # unstated isotope is unconstrained


def test_no_isotope_is_a_leading_zero():
    """`[0C]` demands the absence of a mass number.

    Unlike a real isotope it needs no settled element in the box -- it forbids one span bit -- so
    `[0*]` is legal too.
    """
    assert 'carbon13_methane' not in hits('[0C]')
    assert 'methane' in hits('[0C]')
    assert 'carbon13_methane' not in hits('[0*]')


def test_isotope_needs_an_element():
    """`prim_apply` resolves a mass number against the element's common isotope, so a box with no
    settled element cannot hold one.  The refusal names the position of the digits."""
    with raises(IncorrectSmarts, match='no element symbol'):
        read_smarts('[13]')
    with raises(IncorrectSmarts, match='no element symbol'):
        read_smarts('[13;D2]')


# ----------------------------------------------------------------------------------------------
# COUNTS: D, h, H, x, r, R
# ----------------------------------------------------------------------------------------------
def test_degree_counts_heavy_neighbours():
    assert hits('[O;D0]') == ['water']
    assert hits('[S;D4]') == ['dimethyl_sulfone']


def test_implicit_h_and_total_h_are_different_questions():
    """`h` counts the implicit hydrogens, `H` all of them.  The methyl of `[2H]C` is the case that
    separates them -- three implicit hydrogens, four in total."""
    assert read_smarts('[C;h3]').is_substructure(MOLECULES['deuteriomethane'])
    assert not read_smarts('[C;H3]').is_substructure(MOLECULES['deuteriomethane'])
    assert read_smarts('[C;H4]').is_substructure(MOLECULES['deuteriomethane'])
    assert hits('[C;H0]') == sorted(set(hits('[C;H0]')) - {'methane'}) != []


def test_heteroatom_count():
    assert hits('[C;x2]') == ['acetic_acid', 'benzoic_acid']
    assert 'ethanol' in hits('[C;x1]')


def test_ring_size_and_ring_count():
    assert hits('[C;r3]') == ['cyclopropane']
    assert hits('[C;r5]') == ['cyclopentane', 'pyrrole', 'thf']
    assert 'benzene' in hits('[C;r6]')
    assert 'cyclopropane' in hits('[C;R]') and 'propane' not in hits('[C;R]')
    assert 'propane' in hits('[C;!R]') and 'benzene' not in hits('[C;!R]')


def test_ring_size_zero_is_ring_count_zero():
    """`r0` is `ring_count 0` -- "in no ring" -- which is Daylight's reading of a zero ring size and
    the only one that means anything.  `r1` and `r2` are refused, because a ring of one or two atoms
    does not exist and the box layout has no bit to hold the demand."""
    assert hits('[C;r0]') == hits('[C;!R]')
    for pattern in ['[C;r1]', '[C;r2]']:
        with raises(IncorrectSmarts, match='sizes start at 3'):
            read_smarts(pattern)


# ----------------------------------------------------------------------------------------------
# HYBRIDIZATION, and the one deliberate divergence from V2
# ----------------------------------------------------------------------------------------------
def test_hybridization_is_passed_through_unchanged():
    assert hits('[C;z2]') == ['acetic_acid', 'acetone', 'allene', 'benzoic_acid', 'propene']
    assert hits('[C;z3]') == ['acetonitrile', 'propyne']
    assert 'benzene' in hits('[C;z4]') and hits('[C;z4]') == hits('[C;a]')


def test_z3_is_sp_and_nothing_else():
    """THE DELIBERATE DIVERGENCE.  Measured against chython 2.24: `[C;z3]` there matches allene's
    middle carbon and `[S;z3]` a sulfone's sulfur, because V2 promoted per double bond and capped at
    three.  Here those are `z5` -- two cumulated doubles, no triple -- and `z3` is sp only.

    The reader does not translate.  A template that means "sulfone" and says `z3` is a template to
    re-read, not a string for the lexer to reinterpret.
    """
    assert not read_smarts('[C;z3]').is_substructure(MOLECULES['allene'])
    assert read_smarts('[C;z5]').is_substructure(MOLECULES['allene'])
    assert not read_smarts('[S;z3]').is_substructure(MOLECULES['dimethyl_sulfone'])
    assert read_smarts('[S;z5]').is_substructure(MOLECULES['dimethyl_sulfone'])
    # the nitro group as written here is the charge-separated form, so its nitrogen carries ONE
    # double bond and is `z2` on both sides -- it is the neutral pentavalent spelling that V2 read
    # as `z3` and the core reads as `z5`, and it is here to keep the two spellings from being
    # confused with each other when the 50 `z3` template primitives are ported
    assert read_smarts('[N;+;z2]').is_substructure(MOLECULES['nitrobenzene'])
    assert not read_smarts('[N;+;z3]').is_substructure(MOLECULES['nitrobenzene'])
    assert not read_smarts('[N;+;z5]').is_substructure(MOLECULES['nitrobenzene'])


def test_z5_and_z6_exist_at_all():
    """Six is the storable maximum and seven is where the refusal starts."""
    for value in (1, 2, 3, 4, 5, 6):
        assert read_smarts('[C;z%d]' % value).atom_count_sealed() == 1
    with raises(IncorrectSmarts, match='outside 1..6'):
        read_smarts('[C;z7]')


# ----------------------------------------------------------------------------------------------
# CHARGE, RADICAL, STEREO, MAP NUMBER, MASK
# ----------------------------------------------------------------------------------------------
def test_an_unstated_charge_means_neutral():
    """Both readers agree and it surprises people, so it is pinned: `[N]` does not match a
    quaternary ammonium.  A bracket atom states its charge or states neutrality by omission."""
    ammonium = read_smiles('C[N+](C)(C)C')
    assert not read_smarts('[N]').is_substructure(ammonium)
    assert read_smarts('[N;+]').is_substructure(ammonium)


def test_charge_spellings():
    assert hits('[N;+]') == ['nitrobenzene']
    assert hits('[O;-]') == ['nitrobenzene']
    assert read_smarts('[C;+2]').atom_count_sealed() == read_smarts('[C;++]').atom_count_sealed()
    with raises(IncorrectSmarts, match='outside the storable'):
        read_smarts('[C;+9]')


def test_radical_comes_from_the_extension_tail():
    """`|^1:idx|` after the string, exactly as in CXSMILES, and `idx` is a ZERO-based atom index."""
    q = read_smarts('CC |^1:0|')
    assert q.is_substructure(read_smiles('CC |^1:0|'))
    assert not q.is_substructure(MOLECULES['propane'])   # no radical anywhere


# The charge/radical grid, which the main panel cannot state: it holds one neutral metal, one charged
# molecule and no radical at all.  Local rather than folded into `PANEL` so the primitive tests above
# keep measuring the panel they were written against.
WILDCARD_PANEL = {
    'iron': '[Fe]', 'iron_dication': '[Fe+2]', 'ferrate': '[Fe-]',
    'iron_radical': '[Fe] |^1:0|', 'iron_dication_radical': '[Fe+2] |^1:0|',
    'methane': 'C', 'methyl_cation': '[CH3+]', 'methylene_dication': '[CH2+2]',
    'methyl_radical': '[CH3] |^1:0|', 'methylene_cation_radical': '[CH2+] |^1:0|',
}
WILDCARD_MOLECULES = {k: read_smiles(v) for k, v in WILDCARD_PANEL.items()}


def wildcard_hits(pattern):
    q = read_smarts(pattern)
    return sorted(k for k, m in WILDCARD_MOLECULES.items() if q.is_substructure(m))


def test_any_charge_is_the_withdrawal_of_a_default_not_a_thirteen_way_or():
    """`*` in a bracket frees the CHARGE, and freeing the charge is the whole of it.

    A box is a conjunction of FORBIDDEN bits and the seal neutralises a charge span nobody touched, so
    "any charge" is not a demand for thirteen values -- it is touching the span and forbidding nothing.
    One box, no bits.  Which is also why it composes like every other box bit: a charge stated beside
    it is still forbidding the other twelve, so `[C;*;+2]` is a dication.

    The element half only bites when nothing else settles the element, so `[C;*]` is still carbon.
    `[A]` is unaffected and stays the NEUTRAL any-atom wildcard -- two different wildcards, which is
    what makes both sayable.
    """
    assert wildcard_hits('[M]') == ['iron']                  # unstated charge means neutral, still
    assert wildcard_hits('[M;*]') == ['ferrate', 'iron', 'iron_dication']
    assert wildcard_hits('[C;*]') == ['methane', 'methyl_cation', 'methylene_dication']
    assert wildcard_hits('[C;*;+2]') == ['methylene_dication']
    assert wildcard_hits('[*;+]') == ['methyl_cation']        # a charge beside `*` is not weakened
    assert wildcard_hits('[C,N;*]') == wildcard_hits('[C;*]')  # applies across an OR of elements
    assert wildcard_hits('*') == wildcard_hits('[*]')         # bare and bracketed cannot drift apart
    assert wildcard_hits('[A]') == ['iron', 'methane']
    assert set(wildcard_hits('[A]')) < set(wildcard_hits('[*]'))
    with raises(IncorrectSmarts, match='not a demand'):
        read_smarts('[!*]')                                  # nothing to negate: `*` withdraws


def test_a_radical_is_spellable_in_the_bracket():
    """`^` in a bracket is the field the `|^1:idx|` tail sets, said in place.

    The tail addresses an atom by INDEX, which is fine for a molecule written out once and miserable
    in a query whose atoms ARE the pattern.  No collision with the dative bond, also `^`: a bond token
    is lexed BETWEEN atoms and this one only ever inside a bracket.

    `*` does not free this field.  There are three readings -- radical, not a radical, either -- and a
    token that guessed one of them would make the other two unsayable, so all three are spelled.
    """
    assert wildcard_hits('[C;^]') == ['methyl_radical']       # neutral, as an unstated charge means
    assert wildcard_hits('[C] |^1:0|') == wildcard_hits('[C;^]')      # one field, two spellings
    assert wildcard_hits('[C;!^]') == wildcard_hits('[C]')    # the default, said out loud
    assert wildcard_hits('[C;^,!^]') == ['methane', 'methyl_radical']
    assert wildcard_hits('[M;^]') == ['iron_radical']
    assert wildcard_hits('[M;*;^]') == ['iron_dication_radical', 'iron_radical']
    assert wildcard_hits('[C;*;^]') == ['methyl_radical', 'methylene_cation_radical']
    # the atom token and the bond token are the same character and do not collide
    assert read_smarts('[M;^]^[N;D3]').is_substructure(read_smiles('[Fe]~N(C)(C)C |^1:0|'))


def test_which_atoms_name_no_element_is_a_question_the_boxes_answer():
    """`wildcard_atoms()`: 'any' for an unconstrained element, 'metal' for `[M]`, absent otherwise.

    A rule table needs this to know which atoms of a pattern are shared CONTEXT -- two matches of one
    rule may overlap there and must not overlap on the site being repaired.  The answer is DERIVED from
    the compiled boxes at seal and not from the token that was typed, which is why the answers below
    fall where they do:

    *   `[A]` and `[*]` agree, because they differ in the charge and radical spans and not the element
        one.  Asking "which token was typed" would answer the wrong question -- `*` is the withdrawal
        of the charge default, not an element statement.
    *   `[D2]` is 'any' although nobody typed a wildcard: it constrains no element, and that is the
        whole of what the flag claims.
    *   a list, a negation and `[M,C]` are all absent, because the atom they describe may be a
        carbon.  The metal answer means EVERY box is the 93-metal mask, compared as a mask rather
        than counted, so 93 hand-written elements would not impersonate it.
    """
    assert read_smarts('[A]').wildcard_atoms() == {1: 'any'}
    assert read_smarts('[*]').wildcard_atoms() == read_smarts('[A]').wildcard_atoms()
    assert read_smarts('[M]').wildcard_atoms() == {1: 'metal'}
    assert read_smarts('[D2]').wildcard_atoms() == {1: 'any'}
    assert read_smarts('[A;D2]').wildcard_atoms() == {1: 'any'}
    # the charge and radical spellings are element-silent, so they change nothing
    assert read_smarts('[M;*;^,!^]').wildcard_atoms() == {1: 'metal'}
    for named in ('[C]', '[#6]', '[13C]', '[C,N]', '[!C]', '[M,C]', '[C;*]'):
        assert read_smarts(named).wildcard_atoms() == {}, named
    # every atom is asked, and the stable ids are the query's own
    q = read_smarts('[A:1]-[C:2]~[M:3]')
    assert q.wildcard_atoms() == {1: 'any', 3: 'metal'}


def test_stereo_and_the_unresolved_spelling():
    assert read_smarts('[C@](F)(Cl)Br').atom_count_sealed() == 4
    assert read_smarts('[C@@](F)(Cl)Br').atom_count_sealed() == 4
    log = []
    read_smarts('[C;@?:1](F)(Cl)Br', log)
    assert log and '@?' in log[0]                        # declared, unresolved, dropped, and SAID
    with raises(IncorrectSmarts, match='cannot be negated'):
        read_smarts('[C;!@]')


def test_map_numbers_are_labels():
    q = read_smarts('[C:1][O:2]')
    assert dict(q.map_numbers()) == {1: 1, 2: 2}
    with raises(IncorrectSmarts, match='cannot be negated'):
        read_smarts('[C;!:1]')
    with raises(IncorrectSmarts, match='above the storable'):
        read_smarts('[C:99999]')


def test_mask_is_a_flag_and_not_a_map_number():
    """Masking is a flag of its own, so a masked atom keeps whatever map number the template gave it,
    or none.  A template ported from V2, where masking WAS a reserved map number, must not carry that
    spelling over."""
    q = read_smarts('[C;M:1][O:2]')
    assert sorted(q.masked_atoms()) == [1]
    assert dict(q.map_numbers()) == {1: 1, 2: 2}
    assert sorted(read_smarts('[C;M][O]').masked_atoms()) == [1]
    assert not dict(read_smarts('[C;M][O]').map_numbers())
    with raises(IncorrectSmarts, match='cannot be negated'):
        read_smarts('[C;!M]')


# ----------------------------------------------------------------------------------------------
# THE LOGIC: `;` `,` `&` and juxtaposition
# ----------------------------------------------------------------------------------------------
def test_or_within_a_primitive_type():
    assert hits('[C;D2,D3]') == sorted(set(hits('[C;D2]')) | set(hits('[C;D3]')))
    assert hits('[C;r5,r6]') == sorted(set(hits('[C;r5]')) | set(hits('[C;r6]')))


def test_juxtaposition_binds_tighter_than_or():
    """`[N+,O]` is (nitrogen and +1) or oxygen -- the charge belongs to the alternative that states
    it, and `&` spells the same tight AND explicitly for anyone who wants it written.

    Measured against V2: it splits the body on `;` and then on `,`, builds one element LIST out of the
    pieces and hoists the charge onto the atom, so its `[N+,O]` is `[N;+]` and the oxygen alternative
    is gone.  Precedence is per box here, so a ported `[X+,Y]` means something else than it did.
    """
    assert hits('[N+,O]') == sorted(set(hits('[N;+]')) | set(hits('[O]')))
    assert 'water' in hits('[N+,O]')
    assert hits('[C&D1]') == hits('[C;D1]')
    assert hits('[CD1]') == hits('[C;D1]')


def test_and_low_is_the_documented_separator():
    assert hits('[C;D1;x0;z1]') == hits('[C;z1;x0;D1]')  # order cannot matter


# ----------------------------------------------------------------------------------------------
# BONDS
# ----------------------------------------------------------------------------------------------
def test_implicit_bond_is_single_and_never_aromatic():
    """The dialect's most load-bearing rule and the commonest source of template bugs: an absent
    bond matches order 1 only.  An aromatic bond is written `:`."""
    assert hits('CC') == hits('C-C')
    assert 'benzene' not in hits('CC')
    assert 'benzene' in hits('C:C')
    assert hits('[C;a]:[C;a]') == hits('C:C')


def test_bond_orders():
    assert hits('C=C') == ['allene', 'propene']
    assert hits('C#C') == ['propyne']
    assert hits('C=O') == ['acetic_acid', 'acetone', 'benzoic_acid']


def test_any_bond_matches():
    """RULING.  `~` is the disjunction of every order there is, and `!~` is refused rather than
    compiled into nothing.  V2 lexes `~` to bond order 8, the dative bond, so a template carried over
    verbatim asks a different question here."""
    assert hits('C~C') == sorted(set(hits('C-C')) | set(hits('C=C')) | set(hits('C#C'))
                                 | set(hits('C:C')))
    with raises(IncorrectSmarts, match='forbids every bond'):
        read_smarts('C!~C')


# The dative bond needs molecules `PANEL` does not have, and they must NOT be added to it: the
# panel's expected hit lists are written out in full, so one more member would edit forty of them.
# `~` here is the SMILES spelling of order 8, which is NOT what `~` means in SMARTS -- see below.
DATIVE_PANEL = {
    'borane_ammonia': 'B~N', 'borane_ammonia_misdrawn': 'BN',
    'iron_trimethylamine': '[Fe]~N(C)(C)C', 'trimethylamine': 'CN(C)C',
    'tetramethylammonium': 'C[N+](C)(C)C', 'tetramethylammonium_misdrawn': 'CN(C)(C)C',
    'diborane': 'C[B]1(~[H][B](~[H]1)(C)C)C',
    'propane': 'CCC', 'propene': 'CC=C', 'propyne': 'CC#C', 'benzene': 'c1ccccc1',
    'pyridine': 'c1ccncc1', 'dimethyl_sulfoxide': 'CS(=O)C',
}
DATIVE_MOLECULES = {k: read_smiles(v) for k, v in DATIVE_PANEL.items()}


def dative_hits(pattern):
    q = read_smarts(pattern)
    return sorted(k for k, m in DATIVE_MOLECULES.items() if q.is_substructure(m))


def test_dative_bond_is_spelled_caret():
    """RULING.  `^` is the dative bond, order 8, and the only way to ask for one: `~` is the
    disjunction of every order and the four order tokens are the four covalent ones.  A tree that
    stores coordination as its own order needs a query token for it, because the standardisation
    tables that CREATE order-8 bonds have to be able to find one.

    The two dialects spell it differently on purpose.  In SMILES `~` IS the dative bond -- that is
    how the molecules above are written -- while in SMARTS `~` is "any bond" and `^` is the dative
    one.  Reading `[Fe]~N(C)(C)C` and matching it needs `[M]^[N;D3]`.
    """
    assert dative_hits('[A]^[A]') == ['borane_ammonia', 'diborane', 'iron_trimethylamine']
    assert dative_hits('[B]^[N]') == ['borane_ammonia']
    assert dative_hits('[B]-[N]') == ['borane_ammonia_misdrawn']   # the order tokens exclude it
    assert dative_hits('[M]^[N;D3]') == ['iron_trimethylamine']
    assert dative_hits('[N;D3]^[M]') == ['iron_trimethylamine']    # and the direction is not stated
    assert dative_hits('[H]^[B]') == ['diborane']                  # a three-centre bridge


def test_not_dative_is_expanded_into_the_four_covalent_orders():
    """`!^` cannot be a box the way `!-` is.  "Not order 8" is the set {1, 2, 3, aromatic}, which
    over the order bits is a DISJUNCTION, and a box is a conjunction of forbidden bits -- which is
    why `!~` is refused outright rather than compiled into nothing.  So the reader expands `!^` to
    exactly that disjunction and the caller never has to know: it is `-,=,#,:` and measurably so.
    """
    assert dative_hits('[A]!^[A]') == dative_hits('[A]-,=,#,:[A]')
    assert 'borane_ammonia' not in dative_hits('[A]!^[A]')
    assert 'borane_ammonia_misdrawn' in dative_hits('[A]!^[A]')


def test_any_bond_includes_the_dative_one():
    """`~` is five orders, not four -- so `[B]~[N]` finds the adduct however it was drawn, which is
    what a rule looking for one wants, and `[B]^[N]` is how it says "only the dative drawing"."""
    assert dative_hits('[B]~[N]') == ['borane_ammonia', 'borane_ammonia_misdrawn']
    assert dative_hits('[A]~[A]') == sorted(set(dative_hits('[A]!^[A]'))
                                            | set(dative_hits('[A]^[A]')))


def test_degree_and_heteroatom_count_ignore_a_dative_bond():
    """RULING, and the reason a rule-table row can write `D4` instead of `(-[*])(-[*])(-[*])-[*]`.

    `D` and `x` count SUBSTITUENTS, and a coordination contact is not one.  So `D4` on a nitrogen is
    the claim "four substituents, therefore a formal charge" and cannot be satisfied by a
    three-coordinate donor -- which it could if `D` counted every bond, and `standardize()` would then
    charge a metal-bound amine and grow it a phantom hydrogen.

    `z` counts the same way, so what this test pins is that the three agree.
    """
    assert dative_hits('[N;D4;z1]') == ['tetramethylammonium_misdrawn']
    assert dative_hits('[N;D3;z1;x0]') == ['iron_trimethylamine', 'trimethylamine']
    assert dative_hits('[Fe;D0;x0]') == ['iron_trimethylamine']  # the acceptor's only bond is dative
    # the stored, structural counts are the other answer, and both are right -- `test_features.py`
    # pins that pair.  Here: `heteroatoms_of` sees the iron, `x` does not.
    iron_amine = DATIVE_MOLECULES['iron_trimethylamine']
    nitrogen = next(s for s in iron_amine if iron_amine.element_of(s) == 7)
    assert iron_amine.degree_of(nitrogen) == 4 and iron_amine.heteroatoms_of(nitrogen) == 1


def test_ring_membership_of_a_bond():
    assert hits('C-;@C') == ['cyclohexane', 'cyclopentane', 'cyclopropane', 'thf']
    assert 'cyclohexane' not in hits('C-;!@C')
    assert 'propane' in hits('C-;!@C')


def test_bond_or():
    assert hits('C-,=C') == sorted(set(hits('C-C')) | set(hits('C=C')))
    assert hits('C-,:C') == sorted(set(hits('C-C')) | set(hits('C:C')))


def test_lowercase_atoms_are_aromatic_and_so_is_the_bond_between_them():
    """RULING.  Measured against 2.24: V2's tokeniser collapses the aromatic-atom token onto the plain
    one, so `smarts('c1ccccc1')` there is a single-bonded carbocycle -- it matches CYCLOHEXANE and
    misses benzene.

    Nothing in this codebase writes lowercase SMARTS, so the divergence moves no template.  It does not
    weaken the implicit-bond rule above: that is about atoms written UPPERCASE with an `a` primitive,
    which is how every template spells it.
    """
    assert hits('c1ccccc1') == hits('[C;a]1:[C;a]:[C;a]:[C;a]:[C;a]:[C;a]:1')
    assert 'benzene' in hits('c1ccccc1')
    assert 'cyclohexane' not in hits('c1ccccc1')
    assert hits('cc') == hits('C:C')
    assert hits('[c]') == hits('[C;a]')
    assert hits('cC') == hits('[C;a]-[C]')               # ONE lowercase end is not enough


# ----------------------------------------------------------------------------------------------
# STRUCTURE: branches, ring closures, components
# ----------------------------------------------------------------------------------------------
def test_branches_and_closures():
    assert read_smarts('C(C)(C)C').atom_count_sealed() == 4
    assert read_smarts('C1CCCCC1').bond_count == 6
    assert read_smarts('C%10CCCCC%10').bond_count == 6
    assert hits('C1CC1') == ['cyclopropane']


def test_closure_bond_expression_conflict_is_logged_not_refused():
    """`C-1CC=1` states one bond where the label opens and another where it closes.  Input is
    garbage by default: the reader keeps the opening one, says so, and does not raise."""
    log = []
    read_smarts('C-1CC=1', log)
    assert log and 'the opening one is kept' in log[0]


def test_dot_means_no_bond_stated():
    """As Daylight specifies.  V2 additionally demanded different connected components, so its
    `[O;D1].[O;D1]` missed benzoic acid -- which has two of them on one carbon."""
    q = read_smarts('[O;D1].[O;D1]')
    assert q.bond_count == 0
    assert q.is_substructure(MOLECULES['benzoic_acid'])


def test_a_chain_is_a_monomorphism():
    """`CCC` matches cyclopropane, as Daylight specifies and RDKit agrees.  V2's isomorphism was
    induced -- it rejected any match whose mapped atoms carried a bond the query did not state --
    so it missed all three of these.  The kernel decides this, not the reader: the same query built
    through `QueryContainer`'s journal by hand behaves identically."""
    assert read_smarts('CCC').is_substructure(MOLECULES['cyclopropane'])
    assert read_smarts('CCCC').is_substructure(MOLECULES['cyclopropane']) is False  # only 3 atoms
    assert read_smarts('CCCCC').is_substructure(MOLECULES['cyclopentane'])

    q = QueryContainer()
    a, b, c = q.add_atom(), q.add_atom(), q.add_atom()
    for n in (a, b, c):
        q.atom_primitive(n, 'element', 6)
    q.add_bond(a, b), q.bond_primitive(a, b, 'bond_order', 1)
    q.add_bond(b, c), q.bond_primitive(b, c, 'bond_order', 1)
    assert q.is_substructure(MOLECULES['cyclopropane'])


# ----------------------------------------------------------------------------------------------
# REFUSALS.  Every message ends in a byte offset.
# ----------------------------------------------------------------------------------------------
@mark.parametrize('pattern,fragment', [
    ('', 'no atoms in the string'),
    ('[C;', 'unterminated bracket atom at position 0'),
    ('[]', 'states no primitive'),
    ('[C;Q]', "'Q' names no primitive inside a bracket atom, at position 3"),
    ('C-', 'a bond expression ends the string, at position 1'),
    ('C-=C', 'this query term can never match anything'),
    ('C(.(C))', 'branch opens before any atom at position 3'),
    ('((C))', 'component group opens inside the one at position 0, at position 1'),
    ('()', 'component group at position 0 holds no atom'),
    ('(C', 'component group opens at position 0 and never closes'),
    (')C', 'unbalanced `)` at position 0'),
    ('C(C', 'unbalanced `('),
    ('C1CC', 'ring bond 1 opens at position 1 and never closes'),
    ('C%1C', '`%` needs two digits at position 1'),
    ('C11', 'closes on its own atom'),
    # `#0` is not out of range -- it is the R marker, refused by the seal; see test_r_smirks.py
    ('[C;#119]', 'atomic number 119 is outside 0..118 at position 3'),
    ('[C;D]', '`D` needs a degree at position 3'),
    ('[C;z]', '`z` needs a hybridization at position 3'),
    ('[C;r]', '`r` needs a ring size at position 3'),
    ('[C;x]', '`x` needs a heteroatom count at position 3'),
    ('[C;h]', '`h` needs a hydrogen count at position 3'),
    ('[C:]', 'atom map `:` with no number at position 2'),
    ('-C', 'a bond expression starts a component at position 0'),
    ('C-(C)C', 'bond expression immediately before `(` at position 2'),
    ('C(C-)C', 'bond expression immediately before `)` at position 4'),
    ('C-.C', 'bond expression immediately before `.` at position 2'),
    ('$C', "unexpected '$' at position 0; a query primitive belongs inside a bracket"),
])
def test_refusals(pattern, fragment):
    with raises(IncorrectSmarts) as info:
        read_smarts(pattern)
    assert fragment in str(info.value)


def test_recursive_smarts_is_not_the_dialect():
    """No `$(...)`.  It has never been in chython's SMARTS and the refusal says where."""
    with raises(IncorrectSmarts, match='position'):
        read_smarts('[$(CC)]')


def test_non_ascii():
    with raises(IncorrectSmarts, match='non-ASCII'):
        read_smarts('C—C')


def test_a_term_that_can_never_match_is_refused_at_read_time():
    """The query is sealed before it is returned, so a contradiction is an error where the string
    was written rather than a silent failure to match at the first use.

    Measured against V2: `[C;D1;D2]` there matches propane and benzene, because its body split assigns
    each clause in turn and the last one wins -- a string that says two things reads as the second of
    them with nothing said about the first.
    """
    with raises(IncorrectSmarts):
        read_smarts('[C;!C]')
    with raises(IncorrectSmarts, match='never match'):
        read_smarts('[C;D1;D2]')


def test_the_storable_domain_is_the_limit_and_it_is_stated():
    """9999 is the largest map number an arena or an MDL file can hold, and the refusal names the
    limit at the position where the number was written.  V2 accepts `[C:99999]` and stores the number
    as it stands."""
    with raises(IncorrectSmarts, match='above the storable 9999'):
        read_smarts('[C:99999]')
    assert dict(read_smarts('[C:9999]').map_numbers()) == {1: 9999}


def test_a_leading_zero_is_read_and_not_dropped():
    """The zero is the primitive it looks like, so `[0C]` does not match a 13C.  V2 throws the zero
    away, making its `[0C]` plain carbon."""
    assert 'carbon13_methane' in hits('[C]')
    assert 'carbon13_methane' not in hits('[0C]')


def test_a_branch_needs_an_atom_to_branch_from():
    """`C(.(C))` gives the inner `(` nothing to attach to, so it is refused with the position.  V2
    reads the same string as three carbons, dropping both parentheses.

    A `(` at a COMPONENT position is a component group and not this error, so the refusal is reachable
    only inside a branch."""
    with raises(IncorrectSmarts, match='branch opens before any atom'):
        read_smarts('C(.(C))')


def test_a_lone_parenthesised_fragment_is_one_component_group():
    """`(C)` is a group of one, which matches exactly what `C` matches -- a group constrains fragments
    against each other, and with one fragment there is nothing to constrain."""
    q = read_smarts('(C)')
    assert q.component_groups() == ((frozenset({1}), 0),)
    assert 'methane' in hits('(C)')


# ----------------------------------------------------------------------------------------------
# COMPONENT GROUPING.  `.` says "not bonded" and nothing more; the parentheses say where.
# ----------------------------------------------------------------------------------------------
def test_the_three_grouping_states_are_three_different_questions():
    """The reason grouping had to exist: `.` alone cannot express an intramolecular demand.

    One amine and one alcohol, asked of a molecule carrying both and of a mixture carrying one each.
    Ungrouped matches either; one group matches only the single molecule; two groups only the mixture.
    Without this, a template written for a cyclization would fire across two molecules and one written
    for a coupling would fire inside one.
    """
    intra = read_smiles('NCCCCO')
    inter = read_smiles('NCC.CCO')

    ungrouped = read_smarts('[N;D1].[O;D1]')
    assert ungrouped.is_substructure(intra) and ungrouped.is_substructure(inter)

    together = read_smarts('([N;D1].[O;D1])')
    assert together.is_substructure(intra) and not together.is_substructure(inter)

    apart = read_smarts('([N;D1]).([O;D1])')
    assert not apart.is_substructure(intra) and apart.is_substructure(inter)


def test_a_group_is_reported_per_component_not_per_atom():
    """`set_group` is per atom and the seal reduces it per component -- so a two-atom fragment in a
    group reads back as ONE component carrying it, not two atoms each carrying it."""
    q = read_smarts('(CC).(N)')
    assert q.component_groups() == ((frozenset({1, 2}), 0), (frozenset({3}), 1))


def test_a_branch_inside_a_group_is_still_a_branch():
    """The component-position rule is exactly that: once an atom precedes it, `(` means what it always
    meant, group or no group.  `(CC(C)C.N)` is isobutane and an amine demanded in one molecule."""
    q = read_smarts('(CC(C)C.N)')
    assert q.component_groups() == ((frozenset({1, 2, 3, 4}), 0), (frozenset({5}), 0))
    assert q.is_substructure(read_smiles('NCCC(C)C'))
    assert not q.is_substructure(read_smiles('CC(C)C.N'))


def test_grouping_survives_a_closing_parenthesis_as_a_component_break():
    """`(C)(N)` has no `.` between the fragments and needs none: closing a component group ends the
    component, so what follows cannot bond back into it."""
    q = read_smarts('(C)(N)')
    assert q.component_groups() == ((frozenset({1}), 0), (frozenset({2}), 1))
    assert q.bond_count == 0


# ----------------------------------------------------------------------------------------------
# THE LOG.  Everything dropped says so; nothing is repaired.
# ----------------------------------------------------------------------------------------------
def test_trailing_text_is_reported():
    log = []
    q = read_smarts('C1CC1 and a comment', log)
    assert q.atom_count_sealed() == 3
    assert log and 'not part of it and was ignored' in log[0]


def test_unapplied_extension_fields_are_named():
    log = []
    read_smarts('CC |c:0|', log)
    assert log and 'c:0' in log[0]


def test_radical_field_out_of_range():
    log = []
    read_smarts('CC |^1:5|', log)
    assert log and 'the mark was dropped' in log[0]


def test_unterminated_extension_block():
    log = []
    read_smarts('CC |^1:', log)
    assert log and 'not terminated and was ignored' in log[0]


def test_log_is_optional():
    assert read_smarts('C1CC1 trailing').atom_count_sealed() == 3


# ----------------------------------------------------------------------------------------------
# THE ORACLE.  chython 2.24, in its own interpreter, on the whole table above.
# ----------------------------------------------------------------------------------------------
# The interpreter, the pin, the isolation flag and the identity guards all live in `oracle`, and
# `CHYTHON2_ORACLE` is the one variable that selects the interpreter.

# What V2's lexer refuses outright, measured.  These are not divergences to reconcile: `H` and `[0C]`
# reach box bits V2 cannot spell, `*` and `[c]` are Daylight spellings it rejects, and z5/z6 are the
# hybridizations its validator caps away.
V2_REFUSED = ['[C;H3]', '[C;H0]', '[C;H4]', '[c]', '[*;D4]', '*', '[0*]',
              '[C;z5]', '[C;z6]', '[S;z5]', '[N;z5]', '[C;r1]', '[C;r2]', '[C;+9]', '[C;z7]',
              '[C;#0]', '[C;!@]', '[C;!:1]', '[C;!M]', '[13]', '[13;D2]',
              '[$(CC)]', '[C;Q]', '', '[]', '[C;', 'C-', 'C-=C', ')C', 'C(C', 'C1CC',
              # V2 has no token for the dative bond in either direction: `^` falls through to its
              # SMILES tokeniser ("invalid smiles") and `!^` reaches its bond validator, which refuses
              # it.  So a table that creates order-8 bonds cannot query one there.
              'C^C', 'C!^C', '[B]^[N]',
              # nor a token for the two fields whose default an atom withdraws here: V2's lexer knows
              # no `*` (its element lookup goes looking for a QueryElement named `*`) and no radical
              # primitive, so "a metal of any charge" is unsayable there and "a radical" is sayable
              # only by index, in the tail.
              '[*]', '[*;+]', '[C;*]', '[M;*]', '[!*]', '[C;^]', '[C;!^]', '[M;*;^]', '[C;^,!^]',
              'C%1C', 'C11', '[C;D]', '[C;z]', '[C;r]', '[C;x]', '[C;h]', '[C:]', '-C',
              'C-(C)C', 'C(C-)C', 'C-.C', '$C', '[C;!C]', 'C!~C', 'C—C',
              # V2 has five primitives and one acyclic spelling, no negated element, no `&` and no
              # juxtaposition
              '[C;R]', '[C;r0]', '[!C]', '[C&D1]', '[CD1]']

# What V2 reads DIFFERENTLY, each with its own acceptance test above.  Excluded here because this
# sweep asserts agreement and these are the places where agreement is the wrong answer.
V2_DIVERGES = ['C~C', 'c1ccccc1', 'cc', 'cC', '[C;z3]', '[S;z3]', '[N;z3]',
               '[C;M]', '[C;A]', '[C;M:1][O:2]', '[C;M][O]',
               '[0C]',                                # V2 drops the zero
               '[N+,O]',                              # V2 hoists the charge across the OR
               '[C;D1;D2]', '[C:99999]', '(C)',        # V2 is lenient where the arena is not
               '[O;D1].[O;D1]',                       # V2's component rule
               'CCC', 'CCCC', 'CCCCC', 'C1CC1',       # V2's induced matching
               'CC', 'C-C', 'C-;@C', 'C-;!@C', 'C-,=C', 'C-,:C', 'C=C', 'C#C', 'Br', 'Cl']

_CHILD = r'''
from chython import smarts, smiles

payload = _payload

mols = []
for s in payload['panel']:
    try:
        mols.append(smiles(s))
    except Exception:
        mols.append(None)

out = {}
for p in payload['patterns']:
    try:
        q = smarts(p)
    except Exception as exc:
        out[p] = ['ERR', '%s: %s' % (type(exc).__name__, exc)]
        continue
    row = []
    for m in mols:
        if m is None:
            row.append(None)
            continue
        try:
            # `<=`, not `<`: V2's `__lt__` is a PROPER subgraph and reports False for a query the
            # same size as the molecule, which `is_substructure` here does not.
            row.append(bool(q <= m))
        except Exception:
            row.append(None)
    out[p] = ['OK', row]
_emit(out)
'''


def _ask_oracle(patterns, panel=PANEL):
    """Ask chython 2 which of `patterns` match which of `panel`.

    Resolving the interpreter, pinning the version, passing `-I` and asserting that the child did not
    import the tree under test all live in `chython.core.test.oracle`, which a test enforces is the
    only place the oracle is spawned.  A local copy of those guards with `-I` missing turns a negative
    control into a test that agrees with itself forever.
    """
    return ask(_CHILD, {'panel': list(panel.values()), 'patterns': list(patterns)})


SWEEP = [p for p in ['[C;h2]', '[C;x0]', '[C;a]', '[C;r6]', '[C;z1]', '[C;D2]', '[C;!R]',
                     '[Dy]', '[Ho]', '[Mg]', '[Ag]', '[Ru]', '[Ar]', '[At]', '[Hf]',
                     '[H]', '[M]', '[A]', '[C]', '[#7]', '[N]', '[F,Cl,Br,I]', '[C,N,O;D1;x0]',
                     '[C;D1;x0]', '[N;D1;x0]', '[O;D1;x0]', '[13C]', '[2H]',
                     '[O;D0]', '[S;D4]', '[C;h3]', '[C;x2]', '[C;x1]', '[C;r3]', '[C;r5]',
                     '[C;z4]', '[C;D2,D3]', '[C;r5,r6]', '[N+,O]', '[N;+]', '[O]',
                     '[C;D1]', '[C;D1;x0;z1]', '[C;z1;x0;D1]', '[C;D3]', '[C;+2]', '[C;++]',
                     'Br', 'Cl', 'CC', 'C-C', 'C=C', 'C#C', 'C:C', 'C=O', 'C-;@C', 'C-;!@C',
                     'C-,=C', 'C-,:C', '[C;a]:[C;a]', 'CCC', 'CCCC', 'CCCCC', 'C1CC1',
                     'C(C)(C)C', 'C1CCCCC1', 'C%10CCCCC%10',
                     '[C:1][O:2]', '[C@](F)(Cl)Br', '[C@@](F)(Cl)Br',
                     '[C;a]1:[C;a]:[C;a]:[C;a]:[C;a]:[C;a]:1']
         if p not in set(V2_REFUSED) | set(V2_DIVERGES)]


def _has_induced_embedding(q, mol):
    """Does some embedding of `q` in `mol` add no bond the query did not state?

    V2's matcher only accepts those, so where the core's ONLY embeddings are non-induced V2 says no
    and the pair carries no information about the reader.  This is the narrowest possible exclusion:
    it is per pattern AND per molecule, so `CC` is still compared against every acyclic molecule in
    the panel and only excused on the ring where the difference bites.
    """
    stated = q.bond_count
    for mapping in q.get_mapping(mol):
        atoms = list(mapping.values())
        found = sum(1 for i, n in enumerate(atoms) for m in atoms[i + 1:]
                    if mol.order_of(n, m) is not None)
        if found == stated:
            return True
    return False


def test_oracle_agrees_on_everything_it_can_read():
    """The regression net: 60-odd patterns times the panel, against chython 2.24 in its own
    interpreter.  Skips cleanly when the oracle is not installed, so an unprovisioned checkout still
    runs the suite -- a differential test that cannot be skipped is a differential test nobody runs.

    A disagreement here is not automatically a bug in this reader.  It is a question, and the answer
    is written down in one of the acceptance tests above before the pattern is allowed onto
    `V2_DIVERGES`.  The one exclusion applied INSIDE the loop is the kernel's, not the reader's: a
    pair where the core's every embedding is non-induced tells us nothing, because V2 refuses those
    by construction.  A pattern V2 matches and the core does not is never excused.
    """
    answers = _ask_oracle(SWEEP)
    names = list(PANEL)
    unexpected_refusals, mismatches = [], []
    for pattern in SWEEP:
        status, payload = answers[pattern]
        if status == 'ERR':
            unexpected_refusals.append((pattern, payload))
            continue
        q = read_smarts(pattern)
        for name, v2 in zip(names, payload):
            if v2 is None:
                continue
            mine = q.is_substructure(MOLECULES[name])
            if mine == v2:
                continue
            if mine and not v2 and not _has_induced_embedding(q, MOLECULES[name]):
                continue
            mismatches.append('%s vs %s: V3=%s V2=%s' % (pattern, name, mine, v2))
    assert not unexpected_refusals, ('V2 refused a pattern this sweep believed it could read; move '
                                     'it to V2_REFUSED with a reason: %r' % unexpected_refusals[:5])
    assert not mismatches, mismatches[:20]


def test_z3_is_the_only_z_that_moved():
    """The port table for `z`, measured rather than recalled, over 92 molecules and five elements.

    801 `z` primitives ship in this tree and the epic that translates them needs to know exactly
    which ones are safe to leave alone.  The answer measured here: `z1`, `z2` and `z4` select the
    same atoms in both generations, so every primitive that uses them is safe; `z3` is the only one
    that moved, and it moved in the dangerous direction -- V2's `z3` is a strict SUPERSET, so a
    template carried over verbatim matches strictly less here and fails quietly instead of erroring.

    This is a statement about the two generations' hybridization words, not about the lexer, which
    passes the digit through untouched.  It lives here because it is the sweep that can prove it.
    """
    elements = ('C', 'N', 'S', 'P', 'O')
    agreeing = ['[%s;z%d]' % (e, z) for e in elements for z in (1, 2, 4)]
    answers = _ask_oracle(agreeing + ['[%s;z3]' % e for e in elements], TEMPLATE_PANEL)
    names = list(TEMPLATE_PANEL)

    def both(pattern):
        status, payload = answers[pattern]
        assert status == 'OK', 'V2 could not read %s: %r' % (pattern, payload)
        q = read_smarts(pattern)
        return ({n for n, v in zip(names, payload) if v},
                {n for n in names if q.is_substructure(TEMPLATE_MOLECULES[n])})

    for pattern in agreeing:
        v2, core = both(pattern)
        assert v2 == core, ('%s is documented as meaning the same in both generations but selects '
                            'differently: V2-only %r, core-only %r'
                            % (pattern, sorted(v2 - core), sorted(core - v2)))

    widened = {}
    for element in elements:
        v2, core = both('[%s;z3]' % element)
        assert not core - v2, ('the core matched z3 where V2 did not, so z3 is not a subset after '
                               'all and the port table above is wrong: %r' % sorted(core - v2))
        if v2 - core:
            widened[element] = sorted(v2 - core)
    # named explicitly: a port that only remembers "allene" loses the sulfonyls, which are the bulk
    assert 'S' in widened and 'C' in widened, widened
    assert 'dimethyl_sulfone' in widened['S'] and 'benzenesulfonamide' in widened['S'], widened['S']
    assert 'allene' in widened['C'], widened['C']


def test_oracle_refuses_what_this_file_says_it_refuses():
    """The other half of the same claim: every string on `V2_REFUSED` really is one V2 could not
    read.  Without this the list is a place to hide a disagreement."""
    answers = _ask_oracle(V2_REFUSED)
    read_by_v2 = [p for p, (status, _) in answers.items() if status == 'OK']
    assert not read_by_v2, ('V2 reads these after all, so they belong in the sweep or in an '
                            'acceptance test: %r' % read_by_v2)


# ----------------------------------------------------------------------------------------------
# THE REAL TEMPLATES.  Not spellings chosen to exercise a rule -- the strings this codebase ships.
# ----------------------------------------------------------------------------------------------
_LITERAL = compile_regex(r"smarts\(\s*'([^']{2,300})'")
TEMPLATE_SOURCES = ['algorithms/groups', 'algorithms/standardize', 'algorithms/mapping', 'reactor']


def _template_corpus():
    """Every `smarts('...')` literal shipped in this tree, deduplicated and sorted.

    Read out of the SOURCE rather than out of the rule tables, because the tables hold compiled
    queries and a compiled query can only give its string back through a writer -- which would make
    this a test of the writer.  When the template packages move or go, this returns nothing and the
    tests below skip: an empty corpus is not a passing differential.
    """
    root = Path(__file__).resolve().parent.parent.parent
    out = set()
    for folder in TEMPLATE_SOURCES:
        for path in sorted((root / folder).glob('*.py')) if (root / folder).is_dir() else ():
            out.update(_LITERAL.findall(path.read_text(encoding='utf-8')))
    return sorted(out)


TEMPLATES = _template_corpus()


def test_every_shipped_template_reads():
    """The reader's actual job.  No oracle: a template that does not read is a failure on its own."""
    if not TEMPLATES:
        skip('no `smarts(...)` literals in this tree any more -- point TEMPLATE_SOURCES at the '
             'package that holds the templates now')
    assert len(TEMPLATES) > 100, 'the harvest looks broken, not the reader'
    for pattern in TEMPLATES:
        read_smarts(pattern)


def test_shipped_templates_match_what_v2_matched():
    """The whole shipped corpus against chython 2.24, molecule by molecule.

    This is the test that says the front end can be swapped.  Everything above it pins a rule; this
    one asks the templates themselves whether the rules add up to the same queries.  Same two
    kernel-level exclusions as the sweep: a dot-separated pattern (V2 demanded separate components)
    and a pair with no induced embedding (V2's matcher rejects those).

    What it does NOT cover, stated so nobody reads more into a green tick than is there: only 3 of
    the 47 `z3` templates in the corpus match any panel molecule at all, so agreement on the other
    44 is agreement about nothing.  Their semantics are pinned instead by the port table in
    `test_z3_is_the_only_z_that_moved`, one primitive at a time, which is where the epic that
    translates them should look.
    """
    if not TEMPLATES:
        skip('no `smarts(...)` literals in this tree any more')
    corpus = [p for p in TEMPLATES if '.' not in p]
    answers = _ask_oracle(corpus, TEMPLATE_PANEL)
    names = list(TEMPLATE_PANEL)
    mismatches, unread, agreed_true = [], [], 0
    for pattern in corpus:
        status, payload = answers[pattern]
        if status == 'ERR':
            unread.append((pattern, payload))
            continue
        q = read_smarts(pattern)
        for name, v2 in zip(names, payload):
            if v2 is None:
                continue
            mine = q.is_substructure(TEMPLATE_MOLECULES[name])
            if mine == v2:
                agreed_true += mine
                continue
            if mine and not v2 and not _has_induced_embedding(q, TEMPLATE_MOLECULES[name]):
                continue
            mismatches.append('%s vs %s: V3=%s V2=%s' % (pattern, name, mine, v2))
    assert not unread, ('V2 cannot read a template shipped in this tree, which is a finding about '
                        'the template and not about this reader: %r' % unread[:5])
    assert not mismatches, mismatches[:20]
    # Two readers agree perfectly on a corpus neither one matches, so the count of cells where both
    # said YES is the only number here that measures anything.  It is asserted so that a panel that
    # drifts into irrelevance, or a reader that quietly stops matching, fails instead of passing.
    assert agreed_true > 400, ('only %d template/molecule pairs matched in BOTH readers -- this '
                               'differential is no longer evidence' % agreed_true)
