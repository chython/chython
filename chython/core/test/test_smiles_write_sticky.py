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
"""STICKY SMILES: a string that begins at one named atom and ends at another, so callers can glue text.

The contract is chython 2's, because consumers outside this repository call `mol.sticky_smiles`.  The
traversal is constrained rather than sampled -- `smw_sticky_traverse` blocks the shortest path's later
atoms and defers its successor, and proves termination -- and the tokens are suppressed by the WRITER,
which has the bond's order in hand rather than a finished string to cut characters off.

So the properties below are the specification, and there are only four of them:

1. the written order starts at `left` and ends at `right` (checked through `report=True`, since a
   removed end leaves no token in the text to look at);
2. the text is a whole molecule when nothing is removed, and becomes one again when the caller glues
   the removed atom back on;
3. a kept bond token is never empty and never lies -- `:` where chython 2 wrote `-`;
4. configuration survives detach-and-reattach, and where it CANNOT survive -- the atom loses a
   neighbour for real -- the sign is refused and the atom is reported in `lost` rather than written.

`tries` is accepted and ignored.  Every value frozen here is V3's own: chython 2's traversal is
randomised and unseeded, so it is not a string oracle.  Public compounds only.
"""
from pytest import mark, raises

from chython.core import MoleculeContainer, read_smiles
from chython.core._core import sticky_smiles, write_smiles


# ------------------------------------------------------------------------------------------------
# FIXTURES.  Public compounds and small graph shapes; the shapes are here for the reasons named.
FIXTURES = {
    'ethanol': 'CCO',
    'toluene': 'c1ccccc1C',
    'benzoic_amine': 'OC(=O)c1ccc(N)cc1',       # aromatic ring plus two functional ends
    'tert_butyl_bromide': 'CC(C)(C)Br',         # a degree-4 branch point
    'cyclohexanol': 'C1CCCCC1O',                # a ring with one pendant: the cut-vertex fixture
    'fluoropyridine': 'FC(F)(F)c1ccncc1',
    'propargyl_chloride': 'C#CCCl',             # a triple bond at a terminal
    'alanine': 'N[C@@H](C)C(=O)O',              # tetrahedral centre with a terminal substituent
    'trans_dichloroethene': 'F/C=C/Cl',         # cis/trans with a terminal reference
    'caffeine': 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',  # fused aromatics, four terminals
    'norbornane': 'C1CC2CCC1C2',                # bridged: no atom's removal disconnects it
    'naphthalene': 'c1ccc2ccccc2c1',
    'cyclopropane': 'C1CC1',
    'decyl_bromide': 'CCCCCCCCCCBr',            # a long chain: the path rule has room to be wrong
}


def mol(name):
    return read_smiles(FIXTURES[name])


def severs(m, sid):
    """True when removing `sid` disconnects the molecule -- the precondition `right` must not violate."""
    rest = [x for x in m if x != sid]
    if not rest:
        return False
    seen = {rest[0]}
    stack = [rest[0]]
    while stack:
        for nbr in m.neighbors_of(stack.pop()):
            if nbr != sid and nbr not in seen:
                seen.add(nbr)
                stack.append(nbr)
    return len(seen) != len(rest)


# ------------------------------------------------------------------------------------------------
# 1. WHERE THE STRING STARTS AND ENDS.
@mark.parametrize('name', sorted(FIXTURES))
def test_left_is_the_first_written_atom(name):
    m = mol(name)
    for sid in m:
        text, order, _ = sticky_smiles(m, sid, report=True)
        assert order[0] == sid, (name, sid, text)
        assert len(order) == len(m) and set(order) == set(m)
        assert read_smiles(text).canonical_bytes == m.canonical_bytes


def test_left_is_a_constraint_and_not_a_coincidence():
    # NEGATIVE CONTROL for the test above: the unconstrained writer does NOT start at an arbitrary
    # atom, so `order[0] == sid` is a fact about the sticky walk rather than about every walk.
    m = mol('cyclohexanol')
    default = write_smiles(m, '', True)[1][0]
    assert any(write_smiles(m, '', True)[1][0] != sid for sid in m)
    assert any(sticky_smiles(m, sid, report=True)[1][0] != default for sid in m)


@mark.parametrize('name', sorted(FIXTURES))
def test_right_is_the_last_written_atom(name):
    m = mol(name)
    tried = 0
    for sid in m:
        if severs(m, sid):
            with raises(ValueError, match='cut vertex'):
                sticky_smiles(m, None, sid)
            continue
        tried += 1
        text, order, _ = sticky_smiles(m, None, sid, report=True)
        assert order[-1] == sid, (name, sid, text)
        assert len(order) == len(m) and set(order) == set(m)
        assert read_smiles(text).canonical_bytes == m.canonical_bytes
    assert tried, name       # a fixture where every atom is refused would prove nothing


@mark.parametrize('name', sorted(FIXTURES))
def test_both_ends_at_once(name):
    m = mol(name)
    ids = list(m)
    tried = 0
    for left in ids:
        for right in ids:
            if left == right or severs(m, right):
                continue
            tried += 1
            text, order, _ = sticky_smiles(m, left, right, report=True)
            assert order[0] == left and order[-1] == right, (name, left, right, text)
            assert len(order) == len(m) and set(order) == set(m)
            assert read_smiles(text).canonical_bytes == m.canonical_bytes
    assert tried, name


def test_a_ring_atom_can_be_both_ends():
    # chython 2 requires `right` to be TERMINAL.  The walk does not need that -- what it needs is that
    # nothing hides behind `right` -- so a ring atom is a legal end here, and `...C1` is a good last
    # token.
    m = mol('cyclohexanol')
    text, order, _ = sticky_smiles(m, 1, 3, report=True)
    # the walk goes 1 -> 6 -> 5 -> 4 the long way round, is BLOCKED from 4 to 3, picks up the pendant
    # O, comes back for 2, and reaches 3 last, where the ring closes onto 4.
    assert (text, order) == ('C(C(CC1)O)CC1', (1, 6, 5, 4, 7, 2, 3))
    assert read_smiles(text).canonical_bytes == m.canonical_bytes


def test_a_cut_vertex_right_is_refused_and_says_what_to_pass():
    # cyclohexanol: atom 6 carries the hydroxyl, so O is reachable only through it and no walk can
    # leave it for last.  Cut-vertex-ness, not terminality, is the real precondition.
    m = mol('cyclohexanol')
    assert severs(m, 6)
    with raises(ValueError, match='cut vertex'):
        sticky_smiles(m, 1, 6)
    # NEGATIVE CONTROL: its neighbours in the ring are not cut vertices and are accepted.
    assert not severs(m, 5) and not severs(m, 7)
    assert sticky_smiles(m, 1, 5, report=True)[1][-1] == 5
    assert sticky_smiles(m, 1, 7, report=True)[1][-1] == 7


def test_a_bridged_ring_has_no_cut_vertex_at_all():
    # NEGATIVE CONTROL for `severs`: norbornane is 2-connected, so the refusal above must never fire
    # here, and every ordered pair must be writable.
    m = mol('norbornane')
    assert not any(severs(m, sid) for sid in m)
    for left in m:
        for right in m:
            if left != right:
                assert sticky_smiles(m, left, right, report=True)[1][-1] == right


# ------------------------------------------------------------------------------------------------
# 2. REMOVAL AND GLUE-BACK.  The whole point of the entry point: the text is half a molecule.
SYMBOL = {6: 'C', 7: 'N', 8: 'O', 9: 'F', 17: 'Cl', 35: 'Br', 11: 'Na'}


def terminals(m):
    return [sid for sid in m if len(m.neighbors_of(sid)) == 1]


def token(m, sid):
    """The atom token a caller glues back on.  Organic subset only, which every fixture end is."""
    return SYMBOL[m.atom(sid).element]


@mark.parametrize('name', sorted(FIXTURES))
def test_a_single_bonded_end_glues_back_in_all_four_combinations(name):
    """`keep_bond=True` always glues back.  `keep_bond=False` glues back unless a CONFIGURATION was
    referencing the removed atom, and then the atom is reported in `lost` instead of being lied about."""
    m = mol(name)
    ends = [sid for sid in terminals(m) if m.bond(sid, m.neighbors_of(sid)[0]).order == 1]
    if not ends:
        assert name in ('cyclopropane', 'naphthalene', 'norbornane')    # no terminal atom at all
        return
    for sid in ends:
        head = token(m, sid)
        for keep in (False, True):
            # `sid` as the LEFT end: the caller prepends.
            text, _, lost = sticky_smiles(m, sid, remove_left=True, keep_bond_left=keep, report=True)
            glued = read_smiles(head + text).canonical_bytes == m.canonical_bytes
            assert glued == (keep or not lost), (name, sid, keep, text, lost)
            assert not (keep and lost), (name, sid, text, lost)   # a kept bond loses nothing
            # and as the RIGHT end: the caller appends.  A terminal atom is never a cut vertex, so
            # naming it as `right` is always accepted.
            assert not severs(m, sid)
            text, _, lost = sticky_smiles(m, None, sid, remove_right=True, keep_bond_right=keep,
                                          report=True)
            glued = read_smiles(text + head).canonical_bytes == m.canonical_bytes
            assert glued == (keep or not lost), (name, sid, keep, text, lost)
            assert not (keep and lost), (name, sid, text, lost)   # a kept bond loses nothing


def test_both_ends_are_removed_at_once():
    m = mol('propargyl_chloride')
    text = sticky_smiles(m, 1, 4, remove_left=True, remove_right=True,
                         keep_bond_left=True, keep_bond_right=True)
    assert text == '#CC-'
    assert read_smiles('C' + text + 'Cl').canonical_bytes == m.canonical_bytes
    # NEGATIVE CONTROL: with the tokens kept the same walk writes the whole molecule.
    assert sticky_smiles(m, 1, 4) == 'C#CCCl'


def test_keep_bond_is_ignored_when_the_end_is_not_removed():
    # There is no bond token to keep at an end whose ATOM is still there: the leader has no parent
    # bond, and the last atom's bond to its parent is written by the ordinary path.
    m = mol('ethanol')
    assert sticky_smiles(m, 1, keep_bond_left=True) == sticky_smiles(m, 1) == 'CCO'
    assert sticky_smiles(m, 1, 3, keep_bond_right=True) == sticky_smiles(m, 1, 3) == 'CCO'


def test_a_dropped_bond_drops_its_order_and_the_caller_owns_that():
    # `keep_bond=False` removes the BOND as well as the atom, so a caller who glues a bare atom on
    # gets a SINGLE bond back whatever the original was.  Not a defect and not repairable inside the
    # writer -- the caller asked for the bond to go -- but it is the reason `keep_bond_left=True` is
    # what both in-repo consumers pass, and it is worth pinning.
    m = mol('benzoic_amine')
    dropped = sticky_smiles(m, 3, remove_left=True)                  # 3 is the carbonyl O
    assert dropped == 'C(c1ccc(cc1)N)O'
    assert read_smiles('O' + dropped).canonical_bytes != m.canonical_bytes
    kept = sticky_smiles(m, 3, remove_left=True, keep_bond_left=True)
    assert kept == '=C(c1ccc(cc1)N)O'
    assert read_smiles('O' + kept).canonical_bytes == m.canonical_bytes


def test_the_open_fragment_is_not_a_standalone_molecule():
    # An atom token is computed over the WHOLE molecule, so an open fragment read on its own may get
    # a different implicit-hydrogen count than it had.  That is the contract -- the string is meant to
    # be glued before it is read -- and the test exists so nobody mistakes it for a round trip.
    m = mol('ethanol')
    text = sticky_smiles(m, 1, remove_left=True)
    assert text == 'CO'
    assert read_smiles(text).canonical_bytes == read_smiles('CO').canonical_bytes   # methanol!
    assert read_smiles('C' + text).canonical_bytes == m.canonical_bytes             # ethanol again


# ------------------------------------------------------------------------------------------------
# 3. THE KEPT BOND TOKEN.  Emitted by the writer, which has the bond's order.
@mark.parametrize('name', sorted(FIXTURES))
def test_a_kept_bond_token_is_never_empty(name):
    m = mol(name)
    for sid in terminals(m):
        text = sticky_smiles(m, sid, remove_left=True, keep_bond_left=True)
        assert text[0] in '-=#:~/\\', (name, sid, text)
        # NEGATIVE CONTROL: without `keep_bond_left` the same string starts with an ATOM token.
        assert sticky_smiles(m, sid, remove_left=True)[0] not in '-=#:~/\\'


def test_a_double_and_a_triple_bond_are_written_as_such():
    assert sticky_smiles(mol('propargyl_chloride'), 1,
                         remove_left=True, keep_bond_left=True) == '#CCCl'
    assert sticky_smiles(mol('benzoic_amine'), 3,
                         remove_left=True, keep_bond_left=True) == '=C(c1ccc(cc1)N)O'


def test_an_aromatic_kept_bond_is_written_as_a_colon():
    # THE DIVERGENCE FROM CHYTHON 2, and the reason the token is a writer option.  V2 edits the token
    # list after the fact, and its `_format_bond` returns '' for a single bond and for an aromatic one
    # alike, so `if not smiles[0]: smiles[0] = '-'` writes an aromatic bond as single.
    #
    # The input here is GARBAGE ON PURPOSE -- toluene with its exocyclic bond stored as aromatic --
    # because that is exactly what the project rule says arrives: the record is repaired or reported,
    # never rejected.  The writer reports what the record holds, which is `:`.
    m = mol('toluene')
    m.delete_bond(7, 6)
    m.add_bond(7, 6, 4)
    assert m.bond(7, 6).order == 4
    assert sticky_smiles(m, 7, remove_left=True, keep_bond_left=True) == ':c1ccccc1'
    # NEGATIVE CONTROL: the honest single bond of the untouched molecule writes '-'.
    assert sticky_smiles(mol('toluene'), 7, remove_left=True, keep_bond_left=True) == '-c1ccccc1'


# ------------------------------------------------------------------------------------------------
# 4. STEREO.  Configuration is preserved BY CONSTRUCTION when the bond is kept, and refused -- not
# guessed, not silently dropped -- when the atom really loses a neighbour.
def test_a_tetrahedral_centre_survives_detach_and_reattach():
    # L-alanine, `left` = the methyl carbon on the centre.  The removed atom keeps its slot in the
    # walk and its bond token is written, so the neighbour's written order is unchanged and the sign
    # is recomputed for that order like any other -- note it comes out `@@` where the stored-order
    # string writes `@`, because a sign is a function of the order about to be written and nothing else.
    m = mol('alanine')
    text, order, lost = sticky_smiles(m, 3, remove_left=True, keep_bond_left=True, report=True)
    assert (text, order[0], lost) == ('-[C@@H](C(=O)O)N', 3, ())
    assert read_smiles('C' + text).canonical_bytes == m.canonical_bytes
    # NEGATIVE CONTROL, and the one that matters: the ENANTIOMER is a different molecule, so the
    # equality above is a statement about configuration and not just about connectivity.
    assert read_smiles('N[C@H](C)C(=O)O').canonical_bytes != m.canonical_bytes
    assert read_smiles('C' + text).canonical_bytes != read_smiles('N[C@H](C)C(=O)O').canonical_bytes


def test_a_tetrahedral_sign_is_refused_when_the_bond_goes_too():
    # `remove_left=True, keep_bond_left=False` deletes the atom AND its bond, so the string shows the
    # centre with THREE neighbours.  A parity token there would describe four directions the text does
    # not contain, so the sign is not written and the anchor is reported in `lost`.  This is the same
    # rule as `smw_h_frame_unknown`: an unspellable fact is reported, never approximated.
    m = mol('alanine')
    text, _, lost = sticky_smiles(m, 3, remove_left=True, report=True)
    assert (text, lost) == ('C(C(=O)O)N', (2,))
    assert '@' not in text
    # what the caller gets back by gluing is the FLAT molecule, which is what the text said.
    assert read_smiles('C' + text).canonical_bytes == read_smiles('NC(C)C(=O)O').canonical_bytes
    # NEGATIVE CONTROL: keeping the bond writes the sign and reports nothing.
    kept = sticky_smiles(m, 3, remove_left=True, keep_bond_left=True, report=True)
    assert '@' in kept[0] and kept[2] == ()


def test_a_cis_trans_unit_survives_a_kept_bond_and_is_refused_without_one():
    # The same rule for the other kind of configuration: the removed atom is a REFERENCE of the unit.
    m = mol('trans_dichloroethene')
    text, _, lost = sticky_smiles(m, 1, remove_left=True, keep_bond_left=True, report=True)
    assert (text, lost) == ('/C=C/Cl', ())          # the direction token rides the kept bond
    assert read_smiles('F' + text).canonical_bytes == m.canonical_bytes
    assert read_smiles('F' + text).canonical_bytes != read_smiles(r'F/C=C\Cl').canonical_bytes
    text, _, lost = sticky_smiles(m, 1, remove_left=True, report=True)
    assert (text, lost) == ('C=C/Cl', (2,))
    assert read_smiles('F' + text).canonical_bytes == read_smiles('FC=CCl').canonical_bytes


def test_a_stereocentre_can_be_an_end_itself():
    # `left`/`right` need not be terminal, so the centre itself can be an end.  Nothing is removed
    # here, so nothing may be lost -- the walk is reordered, not the molecule.
    m = mol('alanine')
    for sid in m:
        if sid == 2 or severs(m, sid):
            continue
        text, order, lost = sticky_smiles(m, 2, sid, report=True)
        assert (order[0], order[-1], lost) == (2, sid, ()), (sid, text, lost)
        assert read_smiles(text).canonical_bytes == m.canonical_bytes
    # and the other direction is refused for a reason that has nothing to do with stereo: a branch
    # point in an ACYCLIC molecule is a cut vertex, so no string can end on it.
    assert severs(m, 2)
    with raises(ValueError, match='cut vertex'):
        sticky_smiles(m, 1, 2)


# ------------------------------------------------------------------------------------------------
# 5. REFUSALS.  Every one is "the notation cannot say this", and every message names what to pass.
def test_neither_end_named_is_refused():
    with raises(ValueError, match='either left or right'):
        sticky_smiles(mol('ethanol'))


def test_an_unknown_atom_is_a_key_error():
    with raises(KeyError):
        sticky_smiles(mol('ethanol'), 99)
    with raises(KeyError):
        sticky_smiles(mol('ethanol'), None, 99)


def test_one_atom_cannot_be_both_ends():
    with raises(ValueError, match='same atom'):
        sticky_smiles(mol('ethanol'), 2, 2)


def test_removal_needs_the_end_it_removes():
    with raises(ValueError, match='remove_left needs'):
        sticky_smiles(mol('ethanol'), None, 3, remove_left=True)
    with raises(ValueError, match='remove_right needs'):
        sticky_smiles(mol('ethanol'), 1, remove_right=True)


def test_removing_a_non_terminal_atom_is_refused():
    # `L(A)B` minus `L` is `(A)B`, and an atom carrying a ring closure would leave the digits dangling.
    m = mol('cyclohexanol')
    with raises(ValueError, match='whose degree is 2'):
        sticky_smiles(m, 1, remove_left=True)
    with raises(ValueError, match='whose degree is 3'):
        sticky_smiles(m, 3, 6, remove_right=True)
    # NEGATIVE CONTROL: the one terminal atom of the same molecule is removable.
    assert sticky_smiles(m, 7, remove_left=True, keep_bond_left=True) == '-C1CCCCC1'


def test_removing_both_ends_of_a_two_atom_molecule_is_refused():
    m = read_smiles('CO')
    with raises(ValueError, match='leaves no atom'):
        sticky_smiles(m, 1, 2, remove_left=True, remove_right=True)
    # NEGATIVE CONTROL: removing one end of it is fine.
    assert sticky_smiles(m, 1, 2, remove_left=True, keep_bond_left=True) == '-O'


def test_a_right_end_on_a_salt_is_refused_but_a_left_end_is_not():
    m = read_smiles('CCO.[Na+]')
    assert m.connected_components_count == 2
    with raises(ValueError, match='2 components'):
        sticky_smiles(m, 1, 4)
    text, order, _ = sticky_smiles(m, 1, report=True)
    assert (text, order) == ('CCO.[Na+]', (1, 2, 3, 4))
    assert read_smiles(text).canonical_bytes == m.canonical_bytes


# ------------------------------------------------------------------------------------------------
# 6. THE FROZEN SURFACE.  Consumers outside this repository call the method, not the core function.
def test_the_method_forwards_and_ignores_tries():
    m = mol('ethanol')
    assert m.sticky_smiles(left=1) == sticky_smiles(m, 1) == 'CCO'
    # `tries` is chython 2's retry budget for a randomised traversal.  The walk here is constructed, so
    # every value is accepted and ignored, `tries=0` included.
    assert m.sticky_smiles(left=1, tries=0) == m.sticky_smiles(left=1, tries=10) == 'CCO'
    assert m.sticky_smiles(left=1, right=3, remove_left=True, keep_bond_left=True, tries=1) == '-CO'


def test_the_method_passes_hydrogens_through():
    m = mol('ethanol')
    assert m.sticky_smiles(left=1, hydrogens=True) == '[CH3][CH2][OH]'
    assert m.sticky_smiles(left=1, remove_left=True, keep_bond_left=True,
                           hydrogens=True) == '-[CH2][OH]'
    # NEGATIVE CONTROL: the default writes the organic subset.
    assert m.sticky_smiles(left=1) == 'CCO'


def test_the_result_is_not_canonical_and_does_not_poison_the_cache():
    # The order depends on the atoms the caller named, so this string is not a cache key and must not
    # become one: `str(m)` is the canonical form before and after.  This is the promise the bypass
    # list in `normalize_smiles_spec`'s docstring makes.
    m = mol('cyclohexanol')
    before = str(m)
    assert sticky_smiles(m, 7) != sticky_smiles(m, 1)
    assert str(m) == before == format(m, '')
    assert format(m, 'h') == format(m, 'h')
