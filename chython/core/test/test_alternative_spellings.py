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
"""The four alternative spellings: that each is exactly the value it claims, and that none warns.

A second spelling earns its place by being EXACTLY the thing it forwards to.  One that is nearly right
is worse than none at all, because it ports a consumer silently and wrongly -- so the value tests are
the point of this file.

None of them warns.  `atoms_count` and `len(mol)` are one question asked two ways, so neither is the
one true form, and a library must not print into an application's output for using an API that works.
The silence is asserted rather than assumed: `simplefilter('error')` turns any warning into a failure.

The chython 2 comparison at the foot is an independent witness, not an oracle: it catches a spelling
that disagrees with the library the name comes from.  It runs chython 2 in a subprocess rather than
importing it, so it outlives V2 leaving this tree, and skips when the oracle is not provisioned.

It also records the one place the two libraries disagree -- `[O]`.  That case is excluded from the mass
comparison with the reason written out, not quietly reconciled: a spelling tested against a divergence
would make the divergence the specification.
"""
from math import isclose
from warnings import catch_warnings, simplefilter

from chython.core import read_smiles


# public compounds spanning charge, isotope, implicit hydrogens, radicals and multiple components
COMPOUNDS = ['CC(=O)Oc1ccccc1C(=O)O', 'c1ccccc1O', '[Na+].[Cl-]', 'C[N+](C)(C)C.[Br-]',
             '[13CH4]', 'N#Cc1ccncc1', 'CC(C)(C)OC(=O)N1CCNCC1', 'C', '[O]']

# spelling -> what it must equal, expressed against the other spelling of the same question
EQUIVALENCES = {
    'atoms_count': lambda m: m.atom_count,
    'bonds_count': lambda m: m.bond_count,
    'molecular_charge': int,
    'molecular_mass': float,
}


def test_every_spelling_equals_what_it_forwards_to():
    for smi in COMPOUNDS:
        mol = read_smiles(smi)
        for name, other in EQUIVALENCES.items():
            assert getattr(mol, name) == other(mol), (smi, name)


def test_the_four_spellings_are_alternative_interfaces_and_do_not_warn():
    """Four names chython 2 used that this class answers under a second spelling.

    Not deprecations: `atoms_count` and `len(mol)` are one question asked two ways, so a warning would
    be a library printing into an application's output for using an API that works.
    """
    mol = read_smiles('CC(=O)[O-].[Na+]')
    with catch_warnings():
        simplefilter('error')
        assert mol.atoms_count == len(mol) == 5
        assert mol.bonds_count == mol.bond_count == 3
        assert mol.molecular_charge == int(mol) == 0
        assert mol.molecular_mass == float(mol)


def test_the_two_compatibility_shims_are_gone():
    """`has_atom`/`has_bond` were shims and not spellings -- `has_bond` differed in which exception it
    raised, so a caller had to be edited either way."""
    mol = read_smiles('CCO')
    assert not hasattr(mol, 'has_atom')
    assert not hasattr(mol, 'has_bond')
    assert 1 in mol
    assert mol.order_of(1, 2) is not None


# --- what is deliberately absent -----------------------------------------------------------------

def test_the_computed_names_are_methods_and_not_spellings():
    """`is_radical`, `brutto` and `brutto_formula` are NOT on the list above.

    Each is a first-class property of the container (`test_aggregates.py`), computed rather than
    forwarded.  The distinction the block above draws is between a second spelling and a different
    question, and these are neither -- the same question under the same name.  `element_counts` IS a
    different question and is not a spelling of `brutto`, which is why both exist.
    """
    mol = read_smiles('c1ccccc1O')
    assert mol.is_radical is False
    assert mol.brutto == {'C': 6, 'H': 6, 'O': 1}
    assert mol.brutto_formula == 'C6H6O'
    assert mol.element_counts == {6: 6, 8: 1}


def test_name_is_not_a_spelling_of_title_even_though_it_is_the_same_field():
    """THE HARDEST CASE FOR THE POLICY, AND THE ANSWER IS STILL NO.

    `name` and `title` are the same field -- the record's name line -- and the rename was ruled.  The
    type is not the argument: `title` is `str`, and answers `''` for absent exactly as chython 2's
    `name` did.

    THE SURVIVING REASON IS THE SETTER.  chython 2's `name` is a settable property, so a port carries
    `mol.name = x` as often as it carries `mol.name`.  These containers have no `__dict__` and no
    `name` setter, so a read-only second spelling would make `mol.name = x` raise where the old code
    assigned -- or, on any object that did have a `__dict__`, silently bind a shadow attribute nothing
    reads.  Absent, the port is two one-word edits a grep finds: `mol.name` becomes `mol.title` and
    `mol.name = x` becomes `mol.set_title(x)`.
    """
    mol = read_smiles('CCO')
    mol.set_title(b'ethanol')
    assert mol.title == 'ethanol', 'the surviving spelling, and it is str'
    assert not hasattr(mol, 'name'), \
        'edit the consumer: the spelling is title, and the setter is set_title'


def test_names_whose_meaning_differs_are_not_spellings():
    """The policy, as a test.  `aromatic_rings` exists in chython 2 and is not a rename of `rings`.

    Measured on a molecule where the two answers DIFFER, because the failure is silent.  A consumer
    written against chython 2's `aromatic_rings` that was handed every ring, saturated ones included,
    would draw the wrong molecule and never raise -- and benzene, where both answers are the one same
    ring, would satisfy a forwarding implementation just as happily as a correct one.
    """
    mol = read_smiles('c1ccccc1C1CCCCC1')          # one aromatic ring, one saturated
    assert len(mol.rings) == 2
    assert len(mol.aromatic_rings) == 1, 'a filter over the ring set, not a spelling of rings'


# --- chython 2 as an independent witness ---------------------------------------------------------
#
# Not an oracle. What this catches is a spelling that disagrees with the library the name comes from,
# which is the only thing a second spelling is for. Nothing above depends on it.
#
# Reached through `oracle`: an INSTALLED chython 2 in another interpreter, so this file imports
# no chython 2 and the witness outlives V2 leaving the tree. `from chython import smiles` -- the
# FACADE and not the module path -- is deliberate and stays that way inside the oracle: what is
# being witnessed is the behaviour a consumer of chython 2 actually saw.

V2_VALUES = """
from chython import smiles

out = []
for smi in _payload:
    mol = smiles(smi)
    out.append({'atoms_count': mol.atoms_count, 'bonds_count': mol.bonds_count,
                'molecular_charge': mol.molecular_charge, 'molecular_mass': mol.molecular_mass,
                'numbers': list(mol)})
_emit(out)
"""


def test_the_spellings_agree_with_chython_two():
    from .oracle import ask

    answers = ask(V2_VALUES, COMPOUNDS)
    assert len(answers) == len(COMPOUNDS)
    for smi, old in zip(COMPOUNDS, answers):
        new = read_smiles(smi)
        assert new.atoms_count == old['atoms_count'], smi
        assert new.bonds_count == old['bonds_count'], smi
        assert new.molecular_charge == old['molecular_charge'], smi
        # every number chython 2 numbered an atom is an atom here too, and one that is not an atom
        # in either is not an atom here
        assert all(n in new for n in old['numbers']), smi
        assert max(old['numbers']) + 1 not in new, smi

        if smi == '[O]':
            # NOT COMPARED, AND THE DISAGREEMENT IS chython 2's.  Its SMILES reader re-derives the
            # hydrogen count of a BRACKET atom from valence rules, so `[O]` arrives carrying two
            # implicit hydrogens and masses 18.02 -- where an absent count inside brackets means
            # zero, and this reader gives it zero and masses 15.999.  chython 2 is not even
            # self-consistent about it: `[N]` gets three hydrogens and `[C]` gets none.  The
            # spelling is exact; the molecule being weighed is not the same molecule.
            continue
        # APPROXIMATE, AND ONLY HERE.  That the spelling is exactly `float(mol)` is asserted above by
        # equality; what this line witnesses is that the two libraries weigh the same molecule the
        # same, and they add the masses up in a different grouping -- chython 2 sums atom-plus-its-
        # hydrogens per atom, this sums atoms and hydrogens separately -- so benzene-ol lands
        # 3e-14 apart on the last bits.  Demanding exact equality here would assert an accumulation
        # order neither library promises.
        assert isclose(new.molecular_mass, old['molecular_mass'], rel_tol=1e-12), smi
