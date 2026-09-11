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
"""THE OPERATOR SURFACE, one test per row of chython 2's dunder table.

`MoleculeContainer` answers to twenty-five dunders, and they are the surface a caller reaches for
before any `*_of(n)` accessor.  Each test below is one row of the table, named after the operator, so
a reader can hold the two side by side.

    len(mol)              atom count
    bool(mol)             False only when empty
    iter(mol)             atom NUMBERS, not Atom objects
    x in mol              int -> atom number; str -> element symbol; query/molecule -> substructure
    int(mol)              total formal charge
    float(mol)            molecular mass in daltons
    bytes(mol)            the binary form
    copy.copy(mol)        a copy
    str(mol), format()    canonical SMILES
    repr(mol)             `smiles('...')`, i.e. the same string as an expression that rebuilds it
    mol1 & {n, ...}       the substructure on those numbers
    mol1 - {n, ...}       the substructure without them
    mol1 | mol2           union
    ==, hash              same compound
    <= < >= >             substructure containment
    atom == 'C' / == 6    element symbol or atomic number
    bond == 4             bond order

Three operators chython 2 has are deliberately ABSENT here and asserted absent below: `^`, `~` and
`@` belong to the reactions epic, and a stub that answered them wrongly would be worse than the
`TypeError` a caller gets today.
"""
from copy import copy
from itertools import permutations
from pathlib import Path
from pickle import dumps, loads
from subprocess import run
from sys import executable

from pytest import mark, raises

from chython.core import MoleculeContainer, QueryContainer


def build(atoms, bonds, order=None, **kwargs):
    """`atoms` as symbols or (symbol, kwargs) pairs; returns (molecule, {index: stable id})."""
    m = MoleculeContainer()
    sids = {}
    with m.edit():
        for j in (range(len(atoms)) if order is None else order):
            spec = atoms[j]
            if isinstance(spec, tuple):
                sids[j] = m.add_atom(spec[0], **spec[1])
            else:
                sids[j] = m.add_atom(spec)
        for a, b, o in bonds:
            m.add_bond(sids[a], sids[b], o)
    return m, sids


def chain(*symbols):
    return build(list(symbols), [(i, i + 1, 1) for i in range(len(symbols) - 1)])[0]


PROPANE = ('C', 'C', 'C')
BUTANE = ('C', 'C', 'C', 'C')
ETHANOL = ('C', 'C', 'O')
METHANOL = ('C', 'O')
DIMETHYL_ETHER = ('C', 'O', 'C')

BENZENE = (['C'] * 6, [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 0, 4)])
KEKULE_BENZENE = (['C'] * 6, [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)])
CYCLOHEXANE = (['C'] * 6, [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)])
# hexa-2,4-diene: two cis/trans units in one conjugated chain, which is the fixture the identity
# tests at the bottom of this file need and the reason it is spelled out here.
# THE SP2 CARBONS STATE ONE HYDROGEN AND THE METHYLS STATE NOTHING, which is the smallest honest
# spelling rather than an oversight.  A cis/trans unit is refused when its ANCHOR's hydrogen count is
# unknown -- `-CH=` with two heavy neighbours has three directions or two, and the missing number is
# which -- so every double-bond terminal here has to say.  The methyls are not anchors and nothing
# asks, so they keep `H_UNKNOWN`; giving them their real three hydrogens would be equally true and
# would add two non-stereogenic tetrahedral candidates that no test in this file is about.
_SP2 = ('C', {'implicit_h': 1})
HEXADIENE = (['C', _SP2, _SP2, _SP2, _SP2, 'C'],
             [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1)])
BUTENE = (['C', _SP2, _SP2, 'C'], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])


# ================================================================================================
# len, bool, iter
def test_len_is_the_atom_count():
    """The atom table's length, so implicit hydrogens are not in it."""
    assert len(chain(*PROPANE)) == 3
    assert len(MoleculeContainer()) == 0
    m = chain(*PROPANE)
    assert len(m) == m.atom_count, 'and it is the accessor, not a second count'


def test_bool_is_false_only_when_empty():
    """`if mol:` must not run `len` through a truth table that surprises: an atom is enough."""
    assert not MoleculeContainer()
    assert chain('C')
    assert chain(*PROPANE)


def test_iter_yields_atom_NUMBERS_and_not_atom_objects():
    """The pleasant-looking alternative -- yielding `Atom` views -- would silently change what
    `list(mol)`, `set(mol)` and `dict.fromkeys(mol)` mean, and every loop written as `for n in mol`
    would start handing views to code expecting ints."""
    m, sids = build(list(PROPANE), [(0, 1, 1), (1, 2, 1)])
    assert list(m) == [sids[0], sids[1], sids[2]]
    assert all(isinstance(n, int) for n in m)
    assert set(m) == set(m.atom_numbers)


def test_iter_is_a_snapshot_so_the_loop_survives_an_edit():
    """A view-yielding iterator would go stale mid-loop; a list of numbers does not.  The numbers may
    of course name atoms that no longer exist, which is the caller's problem and a visible one."""
    m, sids = build(list(PROPANE), [(0, 1, 1), (1, 2, 1)])
    seen = []
    for n in m:
        seen.append(n)
        if len(seen) == 1:
            m.delete_atom(sids[2])
    assert seen == [sids[0], sids[1], sids[2]]


# ================================================================================================
# in
def test_in_takes_an_atom_number_as_an_int():
    m, sids = build(list(PROPANE), [(0, 1, 1), (1, 2, 1)])
    assert sids[0] in m
    assert 99 not in m


def test_in_takes_an_element_SYMBOL_as_a_str():
    """And the atomic number spelling is NOT this: `6 in mol` asks about atom number 6."""
    m = chain(*ETHANOL)
    assert 'C' in m
    assert 'O' in m
    assert 'N' not in m
    assert 'Uup' not in m, 'an unknown symbol is a False, not a raise'


def test_in_takes_a_query_or_a_molecule_as_a_substructure_test():
    co = chain(*METHANOL)
    coc = chain(*DIMETHYL_ETHER)
    assert co in coc
    assert coc not in co
    assert co.as_query() in coc
    assert co.as_query() not in chain(*PROPANE)


def test_in_refuses_anything_else():
    with raises(TypeError):
        1.5 in chain(*PROPANE)


# ================================================================================================
# int, float
def test_int_is_the_total_formal_charge():
    """A sum, computed on demand, and deliberately not a stored second truth that an edit could
    leave disagreeing with the atoms."""
    assert int(chain(*PROPANE)) == 0
    m, _ = build([('N', {'charge': 1}), ('O', {'charge': -1}), ('O', {})],
                 [(0, 1, 1), (0, 2, 2)])
    assert int(m) == 0
    m, _ = build([('Na', {'charge': 1})], [])
    assert int(m) == 1


def test_float_is_the_molecular_mass():
    """Ethanol, C2H6O, 46.07 -- so the implicit hydrogens are counted and the tabulated
    abundance-weighted masses are used.  The counts are STATED here rather than derived, because the
    core derives no hydrogens: a mass is only as good as the record it is read from."""
    m, sids = build([('C', {'implicit_h': 3}), ('C', {'implicit_h': 2}), ('O', {'implicit_h': 1})],
                    [(0, 1, 1), (1, 2, 1)])
    assert abs(float(m) - 46.07) < 0.01
    assert float(MoleculeContainer()) == 0.0


def test_float_reads_a_stated_isotope():
    m, sids = build([('C', {'isotope': 13, 'implicit_h': 4})], [])
    assert abs(float(m) - (13.003355 + 4 * 1.00794)) < 0.01
    plain, _ = build([('C', {'implicit_h': 4})], [])
    assert float(m) > float(plain)


# ================================================================================================
# bytes, copy
def test_bytes_is_to_bytes():
    m = chain(*PROPANE)
    assert bytes(m) == m.to_bytes()
    assert MoleculeContainer.from_bytes(bytes(m)).atom_count == 3


def test_copy_copy_is_the_copy_method():
    m = chain(*PROPANE)
    c = copy(m)
    assert c is not m
    assert c.atom_count == 3
    assert list(c.atom_numbers) == list(m.atom_numbers)
    assert c.shares_arena_with(m), 'the arena is immutable, so a copy shares it'


# ================================================================================================
# str, format
def test_str_is_canonical_smiles():
    m = chain(*PROPANE)
    assert str(m) == format(m, '')
    assert str(m) == format(m)


def test_format_forwards_the_spec_and_refuses_an_unknown_key():
    m = chain(*ETHANOL)
    assert format(m, 'i') != '', 'stored slot order is a legal spec'
    with raises(ValueError):
        format(m, 'Q')


def test_format_in_an_f_string():
    m = chain(*PROPANE)
    assert f'{m}' == str(m)
    assert f'{m:i}' == format(m, 'i')


def test_the_smiles_property_is_str():
    """The property spelling of `str(mol)`, and it must not become a SECOND canonicalisation."""
    m = chain(*ETHANOL)
    assert m.smiles == str(m)


def test_repr_is_an_expression_that_rebuilds_the_molecule():
    """`repr` is for pasting back into Python, so it is `smiles('...')` and not the SMILES alone.

    Both halves are asserted, because the useful failure is not "the wrapper is missing" but "the
    string inside it is not this molecule": the call is spelled out and the argument is compared to
    `str`, so a repr that wrapped the STORED-order string would fail here rather than look fine.
    """
    m = chain(*ETHANOL)
    assert repr(m) == "smiles('%s')" % str(m)
    assert repr(m).startswith("smiles('") and repr(m).endswith("')")


def test_repr_of_the_empty_molecule_is_the_CONSTRUCTOR():
    """`smiles('')` is not a molecule, so the one case the general form cannot express says so."""
    assert repr(MoleculeContainer()) == 'MoleculeContainer()'


def test_the_default_string_is_cached_and_the_cache_travels_with_a_copy():
    """`str`, `format('')`, `f'{mol}'` and `.smiles` are ONE cached write, and `copy` inherits it.

    Identity (`is`) and not equality, because the point is that the second read did not recompute: a
    write that produced an equal string every time would pass an `==` check while still costing 6 µs.
    A copy shares the arena, the string is a function of the arena, so recomputing it there would buy
    nothing -- the same argument `canonical_bytes` makes two tests below.
    """
    m = chain(*ETHANOL)
    first = str(m)
    assert str(m) is first
    assert format(m, '') is first
    assert format(m) is first
    assert f'{m}' == first
    assert m.smiles is first
    assert m.copy().smiles is first, 'the string travels with a copy over the same arena'
    # and a non-default spec is NOT served from it, or `i` would answer with the canonical string
    assert format(m, 'i') is not first


MUTATORS = [('set_charge', lambda m, ids: (m.set_charge(ids[2], -1), m.set_hydrogens(ids[2], 0))),
            ('set_hydrogens', lambda m, ids: m.set_hydrogens(ids[2], 0)),
            ('set_isotope', lambda m, ids: m.set_isotope(ids[0], 13)),
            ('set_radical', lambda m, ids: (m.set_radical(ids[2], True),
                                            m.set_hydrogens(ids[2], 0))),
            ('delete_atom', lambda m, ids: m.delete_atom(ids[2])),
            ('add in edit()', lambda m, ids: _grow(m, ids))]


def _grow(m, ids):
    with m.edit():
        n = m.add_atom('N')
        m.add_bond(ids[2], n, 1)


@mark.parametrize('name,mutate', MUTATORS)
def test_every_mutation_invalidates_the_cached_string(name, mutate):
    """The failure mode a cache exists to create: a molecule that changed and a string that did not.

    One case per mutating path reachable from the surface, because the invalidation is not per-path
    code -- it is `_gen`, bumped by `_apply` and by every in-place writer -- and this is the test that
    says so by exhaustion rather than by reading the implementation.  A path added later that forgets
    to bump `_gen` fails here if it is listed, and the list is the point.
    """
    m, ids = build(list(ETHANOL), [(0, 1, 1), (1, 2, 1)])
    for k in ids:
        m.set_hydrogens(ids[k], 3 if k < 2 else 1)
    before = str(m)
    mutate(m, ids)
    after = str(m)
    assert after != before, (name, before, after)
    assert after == format(m, ''), name


def test_a_cached_string_is_NOT_served_inside_an_open_edit_scope():
    """The bug the cache shipped with, for about ten minutes: a stale answer where a refusal is owed.

    `_gen` is bumped when a scope CLOSES, so inside an open one the counter still matches the row
    stored before it opened.  A cache read placed before `_require_clean` therefore hands back the
    pre-scope string for a molecule the caller has just added an atom to -- and that is worse than
    slow, because an uncached `write_smiles` refuses the same call.  Order the two the other way and
    this test fails while every other test in the file still passes, which is exactly why it exists.
    """
    m = chain(*ETHANOL)
    for k, sid in enumerate(m):
        m.set_hydrogens(sid, 3 if k < 2 else 1)
    before = str(m)
    with m.edit():
        m.add_atom('N')
        with raises(RuntimeError):
            str(m)
        with raises(RuntimeError):
            m.smiles
        with raises(RuntimeError):
            format(m, '')
    assert str(m) != before, 'and the scope closing invalidates it'


def test_kekule_and_thiele_invalidate_it_although_the_MOLECULE_is_the_same_compound():
    """The two representation changes, which are exactly the mutations `==` does not see.

    `kekule` and `thiele` leave the compound alone and change how it is spelled, so a cache keyed on
    anything semantic -- `canonical_bytes`, a hash -- would happily serve the wrong string here.  It is
    keyed on `_gen`, which they bump, so the string follows the representation.  Then `thiele` puts it
    back and the ORIGINAL string returns, which is the check that this is invalidation and not just
    change.
    """
    m, ids = build([('C', {'implicit_h': 1})] * 6, [(i, (i + 1) % 6, 4) for i in range(6)])
    aromatic = str(m)
    assert m.kekule().changed
    kekule = str(m)
    assert kekule != aromatic, (aromatic, kekule)
    assert m.thiele().changed
    assert str(m) == aromatic


def test_repr_DOES_NOT_RAISE_on_a_molecule_no_reader_could_answer_for():
    """The property that makes `repr` usable in a debugger, and the only leniency in this class.

    A pending journal makes every read on the container raise -- the arena still holds the pre-scope
    state -- and that is exactly the moment somebody is stepping through `edit()` in a debugger, where
    a raising `repr` replaces the object in the variables pane with a traceback.  So it degrades, and
    the assertion is that the degraded form still NAMES THE REASON: a `repr` that printed a bare
    `<MoleculeContainer>` would send the reader looking for a bug in the wrong place.
    """
    m = MoleculeContainer()
    with m.edit():
        m.add_atom('C')
        text = repr(m)
    assert text.startswith('<MoleculeContainer'), text
    assert 'unwritable' in text and 'RuntimeError' in text, text
    assert 'pending' in text, text
    # and the same molecule, one line later, is perfectly printable
    assert repr(m).startswith('smiles('), repr(m)


# ================================================================================================
# &, -, |
def test_and_is_the_substructure_on_the_given_numbers():
    m, sids = build(list(BUTANE), [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    part = m & {sids[0], sids[1]}
    assert part.atom_count == 2
    assert part.bond_count == 1
    assert list(part.atom_numbers) == [sids[0], sids[1]], 'atom numbers are PRESERVED'


def test_sub_is_the_substructure_without_them():
    m, sids = build(list(BUTANE), [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
    part = m - {sids[3]}
    assert list(part.atom_numbers) == [sids[0], sids[1], sids[2]]
    assert part.bond_count == 2


def test_sub_names_the_numbers_it_does_not_know():
    m = chain(*PROPANE)
    with raises(ValueError) as e:
        m - {99}
    assert '99' in str(e.value)


def test_or_is_the_union_and_remaps_the_right_side():
    m, sids = build(list(PROPANE), [(0, 1, 1), (1, 2, 1)])
    other, _ = build(list(METHANOL), [(0, 1, 1)])
    u = m | other
    assert u.atom_count == 5
    assert u.connected_components_count == 2
    assert list(u.atom_numbers)[:3] == [sids[0], sids[1], sids[2]], 'the left side keeps its numbers'


def test_union_without_remap_refuses_a_collision_and_names_it():
    m = chain(*PROPANE)
    other = chain(*METHANOL)
    with raises(ValueError):
        m.union(other, remap=False)
    other.remap({n: n + 10 for n in other.atom_numbers})
    assert m.union(other, remap=False).atom_count == 5


# ================================================================================================
# ==, hash
def test_eq_ignores_creation_order_and_atom_numbers():
    """`smiles('CCO') == smiles('OCC')`, which is the row Ramil named."""
    cco = chain(*ETHANOL)
    occ, _ = build(['O', 'C', 'C'], [(0, 1, 1), (1, 2, 1)])
    assert cco == occ
    assert not cco != occ
    assert hash(cco) == hash(occ)
    renumbered = cco.copy()
    renumbered.remap({n: n + 50 for n in list(renumbered.atom_numbers)})
    assert cco == renumbered


def test_eq_separates_molecules_the_union_feature_words_cannot():
    """The measurement that closed the question of what `==` may be built on: propane, butane and
    pentane share ONE value of the private union row, because an OR over atoms cannot count them.  An
    `==` built on that row would report the three as one compound."""
    prop, but = chain(*PROPANE), chain(*BUTANE)
    pent = chain('C', 'C', 'C', 'C', 'C')
    assert prop._union_feature_words == but._union_feature_words == pent._union_feature_words
    assert len({prop, but, pent}) == 3
    assert prop != but and but != pent


def test_eq_is_false_and_does_not_raise_against_a_non_molecule():
    m = chain(*PROPANE)
    assert m != 'CCC'
    assert m != 42
    assert m is not None and m != None  # noqa: E711 -- the operator is the thing under test
    assert m != m.as_query()


def test_eq_separates_an_aromatic_ring_from_its_kekule_twin():
    """Two spellings of benzene are two records of two different inputs.  Nothing here normalises a
    representation to make them agree; kekulise both first if that is the question."""
    arom, _ = build(*BENZENE)
    kek, _ = build(*KEKULE_BENZENE)
    assert arom != kek
    assert hash(arom) != hash(kek)
    assert len({arom, kek}) == 2


def test_hash_makes_a_molecule_a_dict_key_and_a_set_member():
    cco = chain(*ETHANOL)
    occ, _ = build(['O', 'C', 'C'], [(0, 1, 1), (1, 2, 1)])
    assert {cco: 'ethanol'}[occ] == 'ethanol'
    assert len({cco, occ, chain(*PROPANE)}) == 2


def test_hash_and_eq_survive_a_bytes_round_trip():
    m, _ = build(*BENZENE)
    back = MoleculeContainer.from_bytes(m.to_bytes())
    assert back == m
    assert hash(back) == hash(m)


def test_hash_and_eq_survive_a_pickle_round_trip():
    """And the cache does not travel in the pickle: `__reduce__` carries the arena's bytes and nothing
    derived, so the canonical form is recomputed in the new process rather than trusted from the old
    one."""
    m, _ = build(*BENZENE)
    back = loads(dumps(m))
    assert back == m
    assert hash(back) == hash(m)


def test_hash_follows_a_mutation_rather_than_going_stale():
    """A molecule is mutable and hashable, which is unusual; the generation counter is what makes it
    sound.  A caller who mutates a molecule already IN a set gets what that always gets in Python: a
    member that can no longer be found."""
    m, sids = build(list(ETHANOL), [(0, 1, 1), (1, 2, 1)])
    before = hash(m)
    m.set_charge(sids[2], -1)
    assert hash(m) != before
    assert m != chain(*ETHANOL)


def test_the_canonical_form_is_cached_and_the_cache_is_not_a_second_truth():
    m, _ = build(*BENZENE)
    first = m.canonical_bytes
    assert m.canonical_bytes is first, 'the same object, so the cache was read'
    c = m.copy()
    assert c.canonical_bytes is first, 'and it travels with a copy over the same arena'


def test_an_empty_molecule_equals_an_empty_molecule():
    assert MoleculeContainer() == MoleculeContainer()
    assert hash(MoleculeContainer()) == hash(MoleculeContainer())
    assert len({MoleculeContainer(), MoleculeContainer()}) == 1


# ================================================================================================
# <= < >= >
def test_the_four_comparisons_are_substructure_containment():
    """`smiles('CO') < smiles('COC')`, the row Ramil named for these."""
    co, coc = chain(*METHANOL), chain(*DIMETHYL_ETHER)
    assert co <= coc
    assert co < coc
    assert coc >= co
    assert coc > co
    assert not coc <= co
    assert not coc < co


def test_the_strict_pair_carries_chython_twos_length_guard():
    """`len` first, so `<` is antisymmetric by construction and the kernel is never asked about a
    fragment that cannot fit."""
    co = chain(*METHANOL)
    same = chain(*METHANOL)
    assert co <= same and same <= co, 'containment holds both ways'
    assert not co < same and not same < co, 'and the strict form holds neither'


def test_containment_is_a_poset_so_two_rings_are_incomparable():
    arom, _ = build(*BENZENE)
    cyc, _ = build(*CYCLOHEXANE)
    assert not arom < cyc
    assert not cyc < arom
    assert not arom <= cyc
    assert not cyc <= arom


def test_a_query_is_accepted_on_the_contained_side_of_all_four():
    co, coc = chain(*METHANOL), chain(*DIMETHYL_ETHER)
    q = co.as_query()
    assert q <= coc
    assert q < coc
    assert coc >= q
    assert coc > q
    assert q <= co, 'the query matches the molecule it was built from'
    assert not q < co, 'but not strictly: same atom count'


def test_a_molecule_is_refused_on_the_pattern_side():
    """`mol <= q` would ask whether a molecule embeds in a pattern, which the kernel cannot answer in
    that direction.  It raises rather than quietly answering the other question."""
    co, coc = chain(*METHANOL), chain(*DIMETHYL_ETHER)
    q = co.as_query()
    with raises(TypeError):
        coc <= q
    with raises(TypeError):
        coc < q
    with raises(TypeError):
        q >= coc
    with raises(TypeError):
        q > coc


def test_the_comparisons_are_not_stereo_aware_and_that_is_documented():
    """`as_query` demands no parity, so one enantiomer contains the other; the alternative -- a query
    carrying a parity built from a molecule -- does not exist yet."""
    left, lsids = build([('C', {'implicit_h': 1}), 'F', 'Cl', 'Br'],
                        [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    right, rsids = build([('C', {'implicit_h': 1}), 'F', 'Cl', 'Br'],
                         [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    left.set_parity(lsids[0], 1)
    right.set_parity(rsids[0], 2)
    assert left <= right and right <= left
    assert left != right, 'while `==` DOES separate them'


# ================================================================================================
# Atom and Bond equality
def test_an_atom_equals_its_symbol_and_its_atomic_number():
    """The idiom is `if atom == 'H'`, and it is what makes `'C' in mol` read as English."""
    m, sids = build(list(ETHANOL), [(0, 1, 1), (1, 2, 1)])
    a = m.atom(sids[0])
    assert a == 'C'
    assert a == 6
    assert int(a) == 6
    assert a != 'O'
    assert a != 8
    assert a != 1.5


def test_an_atom_equals_its_symbol_regardless_of_isotope_and_charge():
    """`atom == 'C'` asks "is this a carbon".  An isotope-aware answer would make it False on a 13-C,
    which no caller means; `atom.isotope` is the question for that."""
    m, sids = build([('C', {'isotope': 13, 'charge': -1, 'radical': True})], [])
    a = m.atom(sids[0])
    assert a == 'C' and a == 6
    assert a.isotope == 13 and a.charge == -1 and a.radical


def test_two_atom_views_compare_by_the_atom_they_name():
    m, sids = build(list(ETHANOL), [(0, 1, 1), (1, 2, 1)])
    assert m.atom(sids[0]) == m.atom(sids[0])
    assert m.atom(sids[0]) != m.atom(sids[1]), 'two carbons, two atoms'
    assert hash(m.atom(sids[0])) == hash(m.atom(sids[0]))
    other = chain(*ETHANOL)
    assert m.atom(sids[0]) != other.atom(list(other.atom_numbers)[0])


def test_a_bond_equals_its_order_and_four_is_the_aromatic_test():
    """An aromaticity test is written `if bond == 4`, and order 4 IS how the arena stores an aromatic
    bond, so this is exact rather than a perception question."""
    m, sids = build(*BENZENE)
    b = m.bond(sids[0], sids[1])
    assert b == 4
    assert int(b) == 4
    assert b != 1
    single = chain(*PROPANE)
    ids = list(single.atom_numbers)
    assert single.bond(ids[0], ids[1]) == 1


def test_two_bond_views_compare_with_their_endpoints_UNORDERED():
    """`bond(1, 2) == bond(2, 1)` returning False would be a trap: a bond has no direction.  The one
    place a pair of atoms IS ordered is `set_wedge`'s narrow/wide, and that is not spelled with a
    Bond."""
    m, sids = build(list(ETHANOL), [(0, 1, 1), (1, 2, 1)])
    assert m.bond(sids[0], sids[1]) == m.bond(sids[1], sids[0])
    assert hash(m.bond(sids[0], sids[1])) == hash(m.bond(sids[1], sids[0]))
    assert m.bond(sids[0], sids[1]) != m.bond(sids[1], sids[2])


# ================================================================================================
# THE OPERATORS THAT ARE NOT HERE, AND THE ONE THAT ARRIVED
def test_xor_and_invert_are_absent_and_meant_to_be():
    """chython 2 spells the CGR `mol1 ^ mol2`.  There is no `CGRContainer` in chython 3 and none is
    planned -- reaction ML consumes `ReactionModelingView` instead -- so this operator is not waiting
    for an epic, it is decided.  A stub answering it from the core would be a wrong answer in place of
    a `TypeError`.

    `~` IS ABSENT FOR ITS OWN REASON: "every single-molecule step" is `mol.react()` with no partner,
    and two spellings of one call is one too many.  `@` is the operator here that does answer -- see
    below."""
    m, other = chain(*PROPANE), chain(*METHANOL)
    with raises(TypeError):
        m ^ other
    with raises(TypeError):
        ~m


def test_matmul_is_a_core_slot_whose_body_arrives_by_injection():
    """`@` is compiled into the core and implemented in `chython.reactions`.

    It has to be BOTH: a special method resolves through the type's slot, and `MoleculeContainer` is
    a `cdef class` that cannot be extended from outside, so the core owns the operator whatever package
    supplies the corpus behind it.  The consequence worth pinning is the failure mode -- in an
    interpreter that never imported `chython.reactions` the operator exists and raises `ImportError`
    naming the package to import, rather than the `TypeError` of an operator that does not exist or the
    empty enumeration of one with no templates.

    A SUBPROCESS WITH THE FACADE STUBBED, and both halves are needed.  A subprocess because
    registration is global and permanent -- any test in this suite that imports `chython.reactions`
    makes these operators work for the rest of the session.  A stub because `import chython.core` runs
    `chython/__init__.py` first, and the facade imports `chython.reactions` for exactly this side
    effect, so the unregistered state cannot be reached with the real facade on the path at all.

    `~` IS ABSENT ALONGSIDE THEM.  "Every single-molecule step" is `mol.react()` with no partner, so
    the operator goes with the distinction it named, along with `oxidize`, `reduce` and `transform`."""
    assert hasattr(MoleculeContainer, '__matmul__')
    for gone in ('__invert__', 'oxidize', 'reduce', 'transform'):
        assert not hasattr(MoleculeContainer, gone), f'{gone} outlived the four-table split'

    script = ("import sys, types\n"
              "stub = types.ModuleType('chython')\n"
              "stub.__path__ = [%r]\n"
              "sys.modules['chython'] = stub\n"
              "from chython.core import read_smiles\n"
              "m = read_smiles('CCO')\n"
              "for f in (lambda: m @ m, lambda: m.react(m), lambda: m.react(),\n"
              "          lambda: m.functional_groups()):\n"
              "    try:\n"
              "        list(f())\n"
              "    except ImportError as e:\n"
              "        assert 'chython.reactions' in str(e), str(e)\n"
              "    else:\n"
              "        raise AssertionError('answered without the corpus registered')\n"
              "assert 'chython.reactions' not in sys.modules\n"
              % str(Path(__file__).resolve().parent.parent.parent))
    result = run([executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# ================================================================================================
# THE IDENTITY, MEASURED.  These are the tests that decide whether `==` may exist at all: a canonical
# form is only an identity if it separates what chemistry separates and merges what chemistry merges,
# and both directions have to be shown.  The oracle for the merges is chython 2's own answer.
def _configure_double_bond(m, x, y, sub_x, sub_y, want):
    """Configure the cis/trans unit on the `x`=`y` bond so that `sub_x` and `sub_y` stand in relation
    `want` (1 or 2).  Returns the stored parity.

    THE ANCHOR IS LOOKED UP, NOT ASSUMED: which end of a double bond carries the unit is a function of
    the arena's slot order, so a fixture that hardcoded one end would pass in the creation order it was
    written in and raise `KeyError` in the next permutation.  The frame is then stated from the
    anchor's end, and reading it from the other
    end swaps `near` with `far`, which a cis/trans relation is invariant under.

    The parity is SEARCHED rather than computed, for `test_smiles_write_stereo.py`'s reason: computing
    it would reimplement the arithmetic under test, so a sign error would cancel.
    """
    units = {u['anchor'] for u in m.stereo_units()}
    if x in units:
        anchor, frame = x, (sub_x, None, sub_y, None)
    else:
        assert y in units, 'no cis/trans unit on this bond'
        anchor, frame = y, (sub_y, None, sub_x, None)
    for parity in (1, 2):
        m.set_parity(anchor, parity)
        if m.translate_stereo(anchor, frame) == want:
            return parity
    raise AssertionError('neither parity gives %r' % (want,))


def _butene(want, order=None):
    m, sids = build(*BUTENE, order=order)
    _configure_double_bond(m, sids[1], sids[2], sids[0], sids[3], want)
    return m


def _hexadiene(first, second, order=None):
    m, sids = build(*HEXADIENE, order=order)
    _configure_double_bond(m, sids[1], sids[2], sids[0], sids[3], first)
    _configure_double_bond(m, sids[3], sids[4], sids[2], sids[5], second)
    return m


def test_cis_and_trans_butene_are_two_compounds():
    cis, trans = _butene(1), _butene(2)
    assert cis != trans
    assert hash(cis) != hash(trans)
    assert len({cis, trans}) == 2


def test_hexa_2_4_diene_has_exactly_THREE_stereoisomers():
    """The oracle is chython 2, which answers three.  FOUR parity combinations exist and only three
    compounds do, because (2E,4Z) and (2Z,4E) are the same molecule read from its two ends.  A test
    asserting "all four differ" would enshrine the opposite bug, so the merge is asserted separately
    below and not left implied by a count."""
    forms = {(a, b): _hexadiene(a, b) for a in (1, 2) for b in (1, 2)}
    assert len(set(forms.values())) == 3
    assert len({hash(m) for m in forms.values()}) == 3


def test_the_two_mixed_hexadienes_are_ONE_compound():
    """(2E,4Z) and (2Z,4E) name the same molecule from opposite ends, so an identity that separated
    them would be reporting one compound as two -- the failure that is easy to mistake for rigour."""
    assert _hexadiene(1, 2) == _hexadiene(2, 1)
    assert hash(_hexadiene(1, 2)) == hash(_hexadiene(2, 1))
    assert _hexadiene(1, 1) != _hexadiene(1, 2)
    assert _hexadiene(2, 2) != _hexadiene(1, 2)


def test_the_identity_does_not_move_with_the_creation_order():
    """The property that makes it an identity rather than a fingerprint of how the record was built.
    Every permutation of six atoms, four configurations, one value each."""
    for want in ((1, 1), (1, 2), (2, 1), (2, 2)):
        seen = {_hexadiene(want[0], want[1], order=list(o))
                for o in permutations(range(6))}
        assert len(seen) == 1, '%r gave %d values over 720 creation orders' % (want, len(seen))


def test_the_identity_survives_a_bytes_round_trip_with_stereo():
    for want in ((1, 1), (1, 2), (2, 2)):
        m = _hexadiene(want[0], want[1])
        assert MoleculeContainer.from_bytes(m.to_bytes()) == m


def test_enhanced_stereo_is_NOT_in_the_identity_yet():
    """The one gap, asserted so it is a known state rather than a surprise: the groups are stored and
    they are not in the canonical form, so a racemate and a single enantiomer of one skeleton compare
    equal.  Delete this test when the canonical form carries them; do not weaken it in place."""
    single, ssids = build([('C', {'implicit_h': 1}), 'F', 'Cl', 'Br'],
                          [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    single.set_parity(ssids[0], 1)
    racemate = single.copy()
    racemate.set_stereo_group(list(racemate.atom_numbers)[0], 1, 1)
    assert racemate.has_stereo_groups and not single.has_stereo_groups
    assert racemate == single, 'a known gap, not a passing grade'
