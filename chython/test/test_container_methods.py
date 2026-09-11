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
"""A pass that acts on one molecule is reachable as a method on it, whichever package holds the body.

The rule this file ratchets: a caller holds a molecule, not a package, so
`chython.chemistry.fix_resonance(mol)` and `mol.fix_resonance()` both work and are ONE body -- the layer
boundary is the library's business and not the caller's.  Tree-wide because it crosses every layer at
once: the bodies live in `core`, `chemistry` and `formats`, and they arrive on a sealed `cdef class` by
injection.

`saturate` is the one deliberate exception and is asserted as one below.
"""

from pytest import raises

from chython import mol, saturate, smiles
from chython.chemistry import calc_implicit, fix_resonance
from chython.core._core import (_set_resonance_fn, _set_sgroup_fns, detached_smiles,
                                molecule_to_inchi, molecule_to_inchikey)
from chython.formats import add_data_sgroup, data_sgroups


#: A four-carbon chain WITH 2D COORDINATES, which the S-group tests need: a `FIELDDISP` anchor is the
#: mean of the referenced atoms' positions, so a molecule read from SMILES has none to give.
BLOCK_2D = """
  chython

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
M  END
"""


# ------------------------------------------------------------------------------------------------
# `chemistry` BODIES REACHED BY INJECTION

def test_fix_resonance_is_a_method_and_the_function_and_they_are_one_body():
    """The pass runs from either spelling, and the log records land on the molecule either way."""
    one, other = smiles('[O-]C=[NH+]C'), smiles('[O-]C=[NH+]C')
    assert other.fix_resonance() is fix_resonance(one) is True
    assert str(one) == str(other)
    assert len(one.log) == len(other.log) == 1


def test_fix_resonance_answers_false_where_no_neutral_form_exists():
    """A nitro group's separated charges are what the drawing meant, so the method leaves them."""
    nitro = smiles('C[N+](=O)[O-]')
    assert not nitro.fix_resonance()
    assert str(nitro) == 'C[N+]([O-])=O'


def test_the_resonance_method_names_the_package_that_registers_it():
    """The refusal a caller sees if the extension is ever reached without `chemistry` imported."""
    molecule = smiles('CCO')
    try:
        _set_resonance_fn(None)
        with raises(ImportError) as e:
            molecule.fix_resonance()
        assert 'chython.chemistry' in str(e.value)
    finally:
        _set_resonance_fn(fix_resonance)
    assert not molecule.fix_resonance()


def test_saturate_is_deliberately_not_a_method():
    """The one pass with no method, and its absence is asserted so a later addition is a decision.

    It is bond perception for a record that gave connectivity and no orders, so its callers are the
    coordinate readers and where those hand their record over is still being designed.
    """
    assert not hasattr(smiles('CCO'), 'saturate')
    assert callable(saturate)


# ------------------------------------------------------------------------------------------------
# `core` BODIES, WHICH NEED NO HOOK -- the method is the body and the function forwards to it, or the
# other way round; either way there is one derivation and not two that can drift.

def test_calc_implicit_is_a_method_and_the_function_forwards_to_it():
    one, other = smiles('CCO'), smiles('CCO')
    n = one.atom_numbers[0]
    assert other.calc_implicit(other.atom_numbers[0]) == calc_implicit(one, n) == 3
    assert one.implicit_h_of(n) == 3


def test_calc_implicit_stores_unknown_rather_than_zero_where_nothing_derives_a_count():
    """`None` is the answer and `H_UNKNOWN` is what gets stored, never a guessed zero -- a caller
    cannot tell an invented 0 from a real one.  This must also never raise: an atom the valence
    collection says nothing about has to survive a repair pass rather than stop it.

    The atom here is pyrrole's nitrogen, the one class left open by design: whether it carries a
    hydrogen is the ring's answer and not a table's, and `kekule()` is what settles it.  A lone metal
    is NOT an example -- `[Fe]` has a free-atom row and answers 0.
    """
    pyrrole = smiles('c1cc[nH]c1')
    n, = (a for a in pyrrole.atom_numbers if pyrrole.element_of(a) == 7)
    assert pyrrole.calc_implicit(n) is None
    assert pyrrole.implicit_h_of(n) is None
    assert smiles('[Fe]').calc_implicit(1) == 0


def test_calc_implicit_recomputes_rather_than_filling_only():
    """The difference from `derive_hydrogens(fill_only=True)`, and the reason to reach for this one
    after an edit: a count already stored is replaced."""
    molecule = smiles('CCO')
    n = molecule.atom_numbers[0]
    molecule.set_hydrogens(n, 0)
    assert molecule.calc_implicit(n) == 3


def test_detached_smiles_is_a_method_and_the_function():
    one, other = smiles('CCOC'), smiles('CCOC')
    cuts = {10: (one.atom_numbers[1], one.atom_numbers[2])}
    ours = other.detached_smiles({10: (other.atom_numbers[1], other.atom_numbers[2])})
    theirs = detached_smiles(one, cuts)
    assert ours.text == theirs.text == 'C%10C'
    assert ours.open_ids == theirs.open_ids == (10,)


def test_detached_smiles_forwards_its_spec_and_its_reserve():
    molecule = smiles('CCOC')
    cuts = {10: (molecule.atom_numbers[1], molecule.atom_numbers[2])}
    assert '[CH3]' in molecule.detached_smiles(cuts, 'h').text
    # a reserved id is withheld from this fragment's own closures, so it cannot collide on a join
    assert 11 not in molecule.detached_smiles(cuts, '', [11]).closure_ids


def test_the_inchi_properties_are_the_functions():
    molecule = smiles('CCO')
    assert molecule.inchi == molecule_to_inchi(molecule) == 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    assert molecule.inchikey == molecule_to_inchikey(molecule) == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


def test_the_inchi_options_stay_on_the_function():
    """`mol.smiles` and `format(mol, spec)` split the same way: the property is the plain answer and
    the keywords live one call out, so a property never grows a parameter."""
    butene = smiles('C/C=C/C')
    assert butene.inchi.endswith('/b4-3+')                            # the double bond's geometry
    assert not molecule_to_inchi(butene, options='-SNon').endswith('/b4-3+')
    assert butene.inchi == molecule_to_inchi(butene)                   # and the property is optionless


# ------------------------------------------------------------------------------------------------
# `formats` BODIES, the one hook that package registers

def test_the_data_label_helpers_are_methods_and_the_functions():
    one, other = mol(BLOCK_2D), mol(BLOCK_2D)
    n = one.atom_numbers[1]
    theirs = add_data_sgroup(one, 'StereoLabel', '(R)', atoms=[n])
    ours = other.add_data_sgroup('StereoLabel', '(R)', atoms=[other.atom_numbers[1]])
    assert ours.field_data == theirs.field_data == '(R)'
    assert ours.disp[:2] == theirs.disp[:2] == (1.0, 0.0)
    assert [r.field_data for r in other.data_sgroups('StereoLabel')] == \
           [r.field_data for r in data_sgroups(one, 'StereoLabel')] == ['(R)']


def test_the_data_label_method_appends_and_survives_both_ctab_versions():
    molecule = mol(BLOCK_2D)
    n, m = molecule.atom_numbers[1], molecule.atom_numbers[2]
    molecule.add_data_sgroup('StereoLabel', '(R)', atoms=[n])
    molecule.add_data_sgroup('NOTE', ['first', 'second'], atoms=[n, m], bonds=[(n, m)])
    assert sorted(r.name for r in molecule.data_sgroups()) == ['NOTE', 'StereoLabel']
    for version in (2000, 3000):
        back = mol(mol(molecule, version=version))
        assert sorted(r.name for r in back.data_sgroups()) == ['NOTE', 'StereoLabel']


def test_the_data_label_methods_name_the_package_that_registers_them():
    molecule = mol(BLOCK_2D)
    try:
        _set_sgroup_fns()
        for call in (lambda: molecule.add_data_sgroup('X', 'y'), lambda: molecule.data_sgroups()):
            with raises(ImportError) as e:
                call()
            assert 'chython.formats' in str(e.value)
    finally:
        _set_sgroup_fns(add_data_sgroup=add_data_sgroup, data_sgroups=data_sgroups)
    assert molecule.data_sgroups() == []
