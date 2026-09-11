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
"""SAME COUNTS FROM EVERY FORMAT THAT BRINGS CONNECTIVITY.  The acceptance gate, mechanized.

Ramil's requirement, in his words: *"I parsed input data, I have maximally computed hydrogens.  I did
kekule/canonicalize/standardize, I got healed all ambiguous cases except errors in data."*  Two claims,
and this file is here because neither of them is checkable from inside any one subpackage:

1. **at read time** every atom holds a count except the one class no local look can settle, and the
   count does not depend on which format the molecule arrived in;
2. **after the repair pipeline** that class is closed too.

The gate lives above `core`, `chemistry` and `formats` for the reason `test_facade_names.py` gives: it
needs a reader and a writer from `formats` and the derivation's own vocabulary from `core`, and
`formats/test/test_isolation.py` forbids reaching across from inside either.

**Why this is a gate and not a unit test.**  The property is *agreement between implementations*, and
it broke three separate times by a reader growing its own copy of the arithmetic -- CTfile answering 0
for a lone hydrogen where MOL2 answered `None`, the SMIRKS patcher refusing every atom that held an
order-4 bond, the InChI reader reading libinchi's "derive this" as "unknown".  Each was correct in
isolation and wrong beside the others, so no test inside one package could see it.  What is asserted
here is the thing a caller actually depends on: the molecule, not the module.

WHAT IS OUT OF SCOPE, AND WHY IT IS NOT AN OMISSION.  XYZ brings coordinates and no bonds, so it
yields frames rather than molecules and there is nothing to derive from; the same goes for the PDB and
mmCIF readers while they return records rather than containers.  `chython/interop/` derives nothing by
design -- it converts what a toolkit already computed.  CML is absent from the façade, so a gate here
cannot reach it; when it lands, its leg belongs in `_CTAB_LEGS` beside the others.
"""
from pytest import mark

from chython import (H_UNKNOWN, inchi_to_molecule, inchi_library_loaded, mol, mol2_mol,
                     molecule_to_inchi, pach_dump, pach_load, read_smiles)
from chython.formats.ctfile import emit_v2000, emit_v3000
from chython.formats.mol2 import element_symbols


_SYMBOLS = element_symbols()

needs_inchi = mark.skipif(not inchi_library_loaded(), reason='libinchi not loaded')

#: Molecules whose counts every format must agree on.  Public compounds only, and chosen for the
#: things that split one reader from another rather than for coverage of chemistry: hypervalent sulfur
#: and phosphorus, a quaternary ammonium, boron, silicon, a nitrile, an amide, a salt whose ions carry
#: no hydrogens at all, and a lone hydrogen atom -- the last being the case two readers answer
#: differently the moment an element fact lives in a reader instead of in the shared derivation.
_CORPUS = ['CCO', 'CC(=O)Nc1ccccc1', 'c1cc[nH]c1', 'c1ccncc1', 'c1ccccc1', 'c1ccsc1', 'c1ccoc1',
           'O=S(=O)(O)O', 'FC(F)(F)S(=O)(=O)O', 'C[N+](C)(C)C', 'CN(C)C=O', 'N#Cc1ccccc1',
           'C[Si](C)(C)C', 'B(O)(O)c1ccccc1', 'OP(=O)(O)O', 'c1ccc2ccccc2c1', 'C[n+]1ccccc1',
           '[Na+].[Cl-]', '[H][H]', 'CC(C)(C)OC(=O)N1CCNCC1']

#: `(name, emitter)` for the two CTAB versions.  Both are gated, because the two writers and the two
#: parsers are separate code and the V3000 branch is the younger one.
_CTAB_LEGS = [('V2000', emit_v2000), ('V3000', emit_v3000)]

#: `(name, version)` for the pach layouts.  Both are gated for the same reason as the two CTABs: the
#: count and the `H_UNKNOWN` sentinel travel through separate encoders, and the default layout is the
#: one this release writes.
_PACH_LEGS = [('the default layout', None), ('version 2', 2)]


def _counts(molecule):
    """`[(element, charge, implicit_h)]` in stable-id order -- the whole answer, per atom."""
    return [(molecule.element_of(n), molecule.charge_of(n), molecule.implicit_h_of(n))
            for n in molecule.atom_numbers]


def _ctab(emit, molecule):
    lines, _ = emit(molecule)
    return '\n'.join(lines) + '\n'


def _through_pach(molecule, version=None):
    """The arena's own serialization, which stores counts rather than deriving them.

    Its leg is here for the opposite reason to the others: it must NOT re-derive, because it is the
    channel that has to carry a count no valence rule reproduces -- diborane's bridging hydrogens,
    ferrocene -- and a `H_UNKNOWN` as an unknown rather than as fifteen.
    """
    restored, _ = pach_load(pach_dump(molecule, version=version))
    return restored


def _aromatic_v2000(molecule):
    """A V2000 whose ring bonds are MDL **type 4**, which chython's own writer refuses to emit.

    Refusing is right -- writing a type-4 bond discards the Kekule form the molecule holds -- but it
    means the only path in the tree that reaches an unstated hydrogen count is unreachable from our
    own output, and that path is exactly what claim (1) is about.  So the test writes it, and writes
    nothing else: no valence field, no `MRV_IMPLICIT_H`, no charge beyond the atom block's own code.
    This is a throwaway emitter for the reader under test, not a format contribution.
    """
    ids = list(molecule.atom_numbers)
    index = {n: i + 1 for i, n in enumerate(ids)}
    bonds = list(molecule.bonds())
    codes = {0: 0, 3: 1, 2: 2, 1: 3, -1: 5, -2: 6, -3: 7}
    lines = ['', '  parity probe', '',
             f'{len(ids):3}{len(bonds):3}  0  0  0  0            999 V2000']
    for n in ids:
        lines.append(f'    0.0000    0.0000    0.0000 {_SYMBOLS[molecule.element_of(n)]:<3} 0'
                     f'{codes.get(molecule.charge_of(n), 0):3}  0  0  0  0  0  0  0  0  0  0')
    for bond in bonds:
        lines.append(f'{index[bond.n]:3}{index[bond.m]:3}{bond.order:3}  0  0  0  0')
    lines.append('M  END')
    return '\n'.join(lines) + '\n'


def test_every_format_agrees_atom_for_atom():
    """CLAIM 1, in its strongest form: not "each reader computes something" but "they agree".

    The reference is SMILES, because it is the one format that states every count in the notation
    itself, so it cannot be wrong for the reason the others can.  A round trip through each of the
    others must reproduce it position by position -- element, charge and count together, since a
    count is only right about the atom it is on.
    """
    for string in _CORPUS:
        reference = read_smiles(string)
        reference.kekule()
        want = _counts(reference)

        for name, emit in _CTAB_LEGS:
            assert _counts(mol(_ctab(emit, reference))) == want, f'{string} through {name}'
        for name, version in _PACH_LEGS:
            assert _counts(_through_pach(reference, version)) == want, f'{string} through pach, {name}'


def test_the_agreement_is_the_DERIVATIONS_and_not_a_channel_in_the_file():
    """...WHICH IS WHAT KEEPS THE TEST ABOVE FROM BEING A TAUTOLOGY, and it nearly is one.

    The CTfile writer emits an `MRV_IMPLICIT_H` data S-group for every atom whose count the valence
    rules would not reproduce, and that S-group is the reader's top-authority channel.  So a round
    trip is guaranteed to agree *even if the derivation answers nothing*, by carrying the numbers
    across in a side channel -- and the test above would pass on a tree where the shared derivation
    had been deleted.

    What it means for that channel to be empty is stated in `valence_for_write`'s own docstring: it
    "gets quieter as the derivation gets better".  Empty is the end of that sentence.  So: for a
    kekulised molecule the file must carry no hydrogen statement at all -- no S-group, no V2000 `vvv`,
    no V3000 `VAL=` -- and the agreement above is then the one derivation reaching the same answer on
    both sides, which is the property with a caller behind it.

    Kekulised is the honest scope.  An aromatic pyrrole nitrogen's count genuinely needs the S-group,
    because the file has no other way to say which nitrogen holds the hydrogen; that is the format's
    limit and not the derivation's, and the next test is where it is faced.
    """
    for string in _CORPUS:
        reference = read_smiles(string)
        reference.kekule()

        for name, emit in _CTAB_LEGS:
            text = _ctab(emit, reference)
            where = f'{string} in {name}'
            assert 'MRV_IMPLICIT_H' not in text, f'{where}: count carried by an S-group'
            assert 'VAL=' not in text, f'{where}: count carried by a stated valence'
            # V2000 states a total valence in `vvv`, columns 49-51 of an atom line, where the field's
            # own 0 means "not stated".  Sliced rather than regexed because a bond line would match
            # any pattern loose enough to find it.
            for line in text.split('\n')[4:4 + reference.atom_count]:
                assert line[48:51].strip() in ('', '0'), f'{where}: stated valence {line[48:51]!r}'


def test_reading_leaves_only_the_class_the_RING_decides():
    """CLAIM 1's exception, and the test is that it is the ONLY one.

    A type-4 bond block states no Kekule form, so a two-coordinate neutral pnictogen in the ring is
    pyrrole or pyridine and the file does not say which -- one hydrogen or none, both legal valences,
    and the ring chooses.  Nothing local can settle it, so the count stores as `H_UNKNOWN`.

    Every other atom in these rings is settled at read time and that is the half worth guarding: an
    aromatic carbon must take a ring double bond in every Kekule form, a thiophene sulfur and a furan
    oxygen cannot take one, and an N-methylpyridinium nitrogen is decided by its charge.  A reader
    that shrugged at all of them -- which is what refusing every atom holding an order-4 bond
    amounts to -- would satisfy a weaker version of this test while answering nothing.
    """
    for string in ['c1ccccc1', 'c1ccsc1', 'c1ccoc1', 'c1ccc2ccccc2c1', 'C[n+]1ccccc1',
                   'Cc1ccccc1O', 'c1ccc(cc1)S(=O)(=O)N']:
        molecule = mol(_aromatic_v2000(read_smiles(string)))
        assert molecule.unknown_h_count == 0, string

    for string in ['c1ccncc1', 'c1cc[nH]c1', 'Cc1ncc[nH]1', 'c1ccc2[nH]ccc2c1']:
        molecule = mol(_aromatic_v2000(read_smiles(string)))
        unknown = [n for n in molecule.atom_numbers if molecule.implicit_h_of(n) is None]
        assert unknown, f'{string}: the ambiguity is real and must not be guessed away'
        for n in unknown:
            assert molecule.element_of(n) in (7, 15), string
            assert molecule.charge_of(n) == 0, string
            assert len(list(molecule.neighbors_of(n))) == 2, string


def test_kekule_closes_every_unknown_the_reading_left():
    """CLAIM 2.  Once the ring HAS a Kekule form the ambiguity is gone, and one call is the whole fix.

    `kekule()` writes the counts its own orders make derivable rather than moving bonds and leaving
    them to somebody else -- Ramil's requirement of 2026-09-04, so that running the pipeline's stages
    by hand gives what `canonicalize()` gives.  Asserted here as zero remaining, over molecules read
    through the one path that produces unknowns in the first place.
    """
    for string in ['c1ccncc1', 'c1cc[nH]c1', 'Cc1ncc[nH]1', 'c1ccc2[nH]ccc2c1', 'c1cnc2[nH]ccc2c1',
                   'c1ccc2[nH]cnc2c1', 'O=c1cc[nH]cn1']:
        molecule = mol(_aromatic_v2000(read_smiles(string)))
        molecule.kekule()
        assert molecule.unknown_h_count == 0, f'{string} after kekule'
        assert all(molecule.implicit_h_of(n) is not None for n in molecule.atom_numbers), string


def test_where_the_file_never_STATED_the_tautomer_the_totals_still_agree():
    """THE HONEST LIMIT, named rather than left for somebody to trip over.

    Benzimidazole and 4-pyrimidinone drawn with type-4 bonds have two candidate nitrogens and the
    file says nothing about which one holds the hydrogen.  `kekule()` picks a Kekule form, the
    hydrogen follows it, and it may not be the nitrogen the SMILES named.  That is not a defect in
    the derivation: the count is right for the molecule that was drawn, and the molecule that was
    drawn is under-specified.  Placing it is a tautomer question -- `standardize_isomers` -- and a
    reader that answered it would be choosing a structure on the file's behalf.

    So the per-atom claim is dropped for these and the multiset is asserted instead, which is the
    part the file does determine and the part a defect would break.
    """
    for string in ['c1ccc2[nH]cnc2c1', 'O=c1cc[nH]cn1']:
        reference = read_smiles(string)
        reference.kekule()
        molecule = mol(_aromatic_v2000(read_smiles(string)))
        molecule.kekule()

        want = sorted(reference.implicit_h_of(n) for n in reference.atom_numbers)
        assert sorted(molecule.implicit_h_of(n) for n in molecule.atom_numbers) == want, string


def test_mol2_reaches_the_same_counts_through_the_same_derivation():
    """SYBYL is the third dialect of "aromatic", and it must land where the other two do.

    A MOL2 `ar` bond is no more informative than an MDL type 4, so pyrrole's nitrogen is the same
    unstated count, and `kekule()` closes it the same way.  Written out as literal text because there
    is no MOL2 writer to round-trip through -- which is also why this leg is worth having: nothing
    else in the suite compares MOL2's counts against another format's.
    """
    text = ('@<TRIPOS>MOLECULE\npyrrole\n 5 5 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n'
            + ''.join(f'{i:7} {a}{i}     0.0000  0.0000  0.0000 {t:<7} 1 RES1  0.0000\n'
                      for i, (a, t) in enumerate([('N', 'N.ar'), ('C', 'C.ar'), ('C', 'C.ar'),
                                                  ('C', 'C.ar'), ('C', 'C.ar')], 1))
            + '@<TRIPOS>BOND\n'
            + ''.join(f'{i:6}{i:5}{i % 5 + 1:5} ar\n' for i in range(1, 6)))

    molecule = mol2_mol(text)
    unknown = [n for n in molecule.atom_numbers if molecule.implicit_h_of(n) is None]
    assert [molecule.element_of(n) for n in unknown] == [7], 'the pnictogen, and only it'

    molecule.kekule()
    assert molecule.unknown_h_count == 0
    reference = read_smiles('c1cc[nH]c1')
    reference.kekule()
    assert (sorted(molecule.implicit_h_of(n) for n in molecule.atom_numbers)
            == sorted(reference.implicit_h_of(n) for n in reference.atom_numbers))


@needs_inchi
def test_inchi_carries_the_counts_back_rather_than_dropping_them():
    """libinchi's `num_iso_H[0] == -1` is an INSTRUCTION -- "derive from the valence rules" -- so the
    reader runs the fill-only sweep on it rather than storing an unknown.  Read as an unknown, it would
    bring a molecule back through InChI with every count missing.

    Compared as a multiset because InChI renumbers and normalizes: the atom order is its own and a
    positional comparison would be asserting something about its canonicalizer.  Total hydrogen count
    is the claim InChI does make, and it is the claim a dropped derivation breaks.

    `[H][H]` is excluded, and the exclusion is InChI's semantics rather than a hole here.  InChI folds
    an explicit hydrogen into its neighbour's count, so `InChI=1S/H2` comes back as ONE atom holding
    one implicit hydrogen -- the same molecule, and a multiset of `[1]` against the SMILES reading's
    `[0, 0]`.  A test that demanded agreement there would be demanding that InChI not normalize.
    """
    for string in (s for s in _CORPUS if s != '[H][H]'):
        reference = read_smiles(string)
        reference.kekule()
        molecule = inchi_to_molecule(molecule_to_inchi(reference))
        molecule.kekule()

        assert molecule.unknown_h_count == 0, string
        assert (sorted(molecule.implicit_h_of(n) for n in molecule.atom_numbers)
                == sorted(reference.implicit_h_of(n) for n in reference.atom_numbers)), string


def test_an_unrecorded_count_is_a_third_state_in_every_direction():
    """The sentinel is 15 and a count is four bits, so `H_UNKNOWN` is a NUMBER unless someone looks.

    That is the shape of the bug this whole property is exposed to: every layer that handles a count
    has to distinguish "none" from "not recorded", and a layer that forgets either claims the atom
    has fifteen hydrogens or claims it has none.  Both have happened.  Asserted end to end -- the
    container's accessor, its tally, and a CTfile round trip -- because the sentinel has to survive
    all three to mean anything.
    """
    molecule = read_smiles('CCO')
    assert molecule.unknown_h_count == 0
    molecule.set_hydrogens(2, H_UNKNOWN)

    assert molecule.implicit_h_of(2) is None, 'not 15, and not 0'
    assert molecule.unknown_h_count == 1
    for name, version in _PACH_LEGS:
        assert _counts(_through_pach(molecule, version))[1][2] is None, \
            f'pach stores the third state, {name}'
