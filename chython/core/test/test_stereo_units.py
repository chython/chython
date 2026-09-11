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
import re

import pytest

from chython.core import MoleculeContainer, _core

_SYMBOL = re.compile(r'[A-Z][a-z]?')


def _symbols(atoms):
    """'CFClBr' -> ['C', 'F', 'Cl', 'Br']. Never iterate the string: 'Cl' is two characters."""
    return _SYMBOL.findall(atoms)


def _mol(*, atoms, bonds, isotopes=None, charges=None, hydrogens=None):
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e,
                           isotope=0 if isotopes is None else isotopes[k],
                           charge=0 if charges is None else charges[k],
                           implicit_h=0 if hydrogens is None else hydrogens[k])
                for k, e in enumerate(_symbols(atoms))]
        for i, j, o in bonds:
            m.add_bond(sids[i], sids[j], o)
    return m, sids


# The core never derives an implicit hydrogen count -- `add_atom(implicit_h=...)` is the only
# source of one (test_derive.py:38) -- so every record below that needs a hydrogen direction
# states it. A carbon written with three heavy neighbours and no `hydrogens` has THREE
# directions here, not four, which is why the dichloromethane record is refused twice over.
#
# `_mol` PASSES `implicit_h=0` AND NOT `None` FOR AN OMITTED `hydrogens`, which is a statement and
# not a default.  `add_atom`'s own default is now `H_UNKNOWN`, and perception refuses an anchor whose
# count is unknown -- correctly, since three heavy neighbours plus an unknown hydrogen is a
# stereocentre or is not and the missing number is precisely which.  Every fixture here that omits
# `hydrogens` is asserting on a DIRECTION COUNT, so it needs the number to exist; the zero is what
# the assertions were computed against and it is now said out loud.  A fixture that wanted a real
# hydrogen passes `hydrogens=` -- that list is the only thing that ever meant a count.


def test_bromochlorofluoromethane_is_a_tetra_candidate():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 0 and u['anchor'] == sids[0] and u['n_refs'] == 4
    # three sigma neighbours in CSR ascending order, then the implicit H last
    assert u['refs'] == (sids[1], sids[2], sids[3], None)


def test_dichloromethane_is_not_a_candidate():
    m, sids = _mol(atoms='CClCl', bonds=[(0, 1, 1), (0, 2, 1)])
    assert m.stereo_units() == []


def test_dichloromethane_with_its_hydrogens_stated_is_a_candidate():
    """Four directions is the whole of this task's rule, and CH2Cl2 has four.

    The pair of chlorines is interchangeable, so the unit is not stereogenic -- but that is
    decided by the automorphism group, not here. The local test refuses only what it
    is certain about (see the fragment comment), and two sigma neighbours are not that case:
    `[2H]C([H])(Cl)Br` differs from CH2Cl2 only in what hangs off the repeated direction.
    """
    m, sids = _mol(atoms='CClCl', hydrogens=[2, 0, 0], bonds=[(0, 1, 1), (0, 2, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['refs'] == (sids[1], sids[2], None, None)


def test_deuterium_makes_a_candidate():
    # [2H]C([H])(Cl)Br -- chiral by isotope alone; V2 returns [] here
    m, sids = _mol(atoms='CHHClBr', isotopes=[0, 2, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assert len(m.stereo_units()) == 1


def test_both_hydrogens_of_a_deuterated_centre_are_named():
    """[2H]C([H])(Cl)Br -- the case Ruling F26 exists for.

    Its only distinguishing feature is which of the two hydrogens is which, so a rule that
    erased hydrogen identity could say `candidate` here but could never say WHICH configuration:
    two of the four directions would be the same `None`. Naming them makes the record
    expressible, and the two hydrogen directions are ordered between themselves by ascending
    slot, exactly as the heavy ones are.

    Heavy directions first, whatever the slot order: the chlorine and bromine were declared
    after both hydrogens, and they still come first.
    """
    m, sids = _mol(atoms='CHHClBr', isotopes=[0, 2, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    # (Cl, Br, 2H, 1H): heavy slots ascending, then hydrogen slots ascending
    assert units[0]['refs'] == (sids[3], sids[4], sids[1], sids[2])
    assert units[0]['unnamed_mask'] == 0


def test_an_explicit_hydrogen_keeps_the_position_an_implicit_one_had_and_gains_a_name():
    """CHFClBr twice: hydrogen implied, then hydrogen drawn. Same ORDER, and now a name.

    This is the property that decides where a hydrogen direction sorts. Explicitness is a
    drawing choice -- parsers, standardisation and the depiction layer add and drop explicit
    hydrogens -- so a parity stored against this list has to keep meaning the same
    configuration across that change. Fixing the hydrogen direction's POSITION (after every
    heavy slot) is what buys that; erasing its identity is not needed for it, and costs the
    deuterated case above. So: the three heavy directions sit in the same three positions in
    both records, and the fourth is `None` when the hydrogen is implied and the hydrogen's own
    stable id when it is drawn.

    A rule that sorted the hydrogen by its own slot -- the naive reading -- would put it
    anywhere in the list, since the H's index depends on when it was added, and every
    add/remove-explicit-H would silently re-base every stored parity in the molecule.
    """
    implied, isids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                          bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    drawn, dsids = _mol(atoms='CFClBrH',
                        bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    iu = implied.stereo_units()
    du = drawn.stereo_units()
    assert len(iu) == 1 and len(du) == 1
    assert iu[0]['refs'] == (isids[1], isids[2], isids[3], None)
    assert du[0]['refs'] == (dsids[1], dsids[2], dsids[3], dsids[4])
    assert iu[0]['n_refs'] == du[0]['n_refs'] == 4


def _popcount(mask):
    return bin(mask).count('1')


def test_the_unnamed_direction_slots_are_recorded():
    """`unnamed_mask` says WHICH slots hold a direction with no atom: implicit H and the lone pair.

    The automorphism filter rejects a direction list with two or more of them -- an implicit
    hydrogen is always protium and a lone pair is unique, so two unnamed directions are necessarily
    indistinguishable. That rule is exceptionless only because an explicit hydrogen is NAMED
    (Ruling F26); count it here and `[2H]C([H])(Cl)Br` would report 2 and be thrown away.

    Ruling F41 made this a 4-bit mask rather than a count, for the bond kinds below. For an ATOM
    kind the two spellings carry the same information, because the unnamed directions are the
    tail of the list -- so both halves of that guarantee are asserted here: the exact mask, and
    its popcount.
    """
    # toluene: the methyl carbon is 1 heavy direction and 3 implicit hydrogens
    toluene, tsids = _mol(atoms='CCCCCCC', hydrogens=[3, 0, 1, 1, 1, 1, 1],
                          bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2),
                                 (4, 5, 1), (5, 6, 2), (6, 1, 1)])
    units = toluene.stereo_units()
    assert len(units) == 1 and units[0]['anchor'] == tsids[0]
    assert units[0]['unnamed_mask'] == 0b1110 and _popcount(units[0]['unnamed_mask']) == 3

    implied, _ = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                      bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert implied.stereo_units()[0]['unnamed_mask'] == 0b1000
    assert _popcount(implied.stereo_units()[0]['unnamed_mask']) == 1

    drawn, _ = _mol(atoms='CFClBrH',
                    bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assert drawn.stereo_units()[0]['unnamed_mask'] == 0
    assert _popcount(drawn.stereo_units()[0]['unnamed_mask']) == 0

    # the sulfur lone pair is unnamed too: dimethyl sulfoxide is 2 sigma + 1 pi + 1 pair
    dmso, dsids = _mol(atoms='SOCC', hydrogens=[0, 0, 3, 3],
                       bonds=[(0, 1, 2), (0, 2, 1), (0, 3, 1)])
    assert dmso.unit_of(dsids[0])['unnamed_mask'] == 0b1000
    assert _popcount(dmso.unit_of(dsids[0])['unnamed_mask']) == 1


def test_phosphine_oxide_is_a_candidate():
    # a pi neighbour contributes one direction, not two
    m, sids = _mol(atoms='POCCC', bonds=[(0, 1, 2), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['anchor'] == sids[0]


def test_ammonium_is_a_candidate():
    m, sids = _mol(atoms='NCCCC', charges=[1, 0, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assert len(m.stereo_units()) == 1


def test_sulfoxide_lone_pair_counts():
    # 2 sigma + 1 pi + lone pair == 4, on sulfur only
    m, sids = _mol(atoms='SOCC', bonds=[(0, 1, 2), (0, 2, 1), (0, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['refs'][3] is None


def test_the_ylide_form_of_a_sulfoxide_is_the_same_candidate():
    """[O-][S+](C)C -- dimethyl sulfoxide drawn as an ylide instead of with S=O.

    The lone pair is counted from the sulfur's electron budget, not from a bond-order pattern,
    so both drawings leave one pair and both reach four directions. An implementation that
    keyed the lone pair off `order == 2` would see three directions here and lose the centre.
    """
    m, sids = _mol(atoms='SOCC', charges=[1, -1, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['refs'] == (sids[1], sids[2], sids[3], None)


def test_sulfonium_is_a_candidate():
    # trimethylsulfonium: 3 sigma + one lone pair from the remaining electron pair
    m, sids = _mol(atoms='SCCC', charges=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert len(m.stereo_units()) == 1


def test_dimethyl_sulfide_is_not_a_candidate():
    """Sulfur here has two lone pairs, and they contribute ONE direction between them.

    Counting both would give 2 sigma + 2 lone pairs == 4 and invent a centre. Two lone pairs on
    one atom are the same direction twice over -- nothing can ever tell them apart -- so a
    record that needs both to reach four is not a candidate.
    """
    m, sids = _mol(atoms='SCC', bonds=[(0, 1, 1), (0, 2, 1)])
    assert m.stereo_units() == []


def test_trimethylamine_lone_pair_does_not_count():
    m, sids = _mol(atoms='NCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.stereo_units() == []


def test_trimethylphosphine_lone_pair_does_not_count():
    m, sids = _mol(atoms='PCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.stereo_units() == []


def test_symmetric_sulfone_is_not_a_candidate():
    m, sids = _mol(atoms='SOOCC', bonds=[(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)])
    assert m.stereo_units() == []


def test_isotopic_sulfone_is_a_candidate():
    # the two oxygens stopped being identical, so no element whitelist can express this
    m, sids = _mol(atoms='SOOCC', isotopes=[0, 16, 18, 0, 0],
                   bonds=[(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)])
    assert len(m.stereo_units()) == 1


def test_sulfonimidoyl_is_a_candidate():
    m, sids = _mol(atoms='SONCC', bonds=[(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)])
    assert len(m.stereo_units()) == 1


def test_two_pi_neighbours_that_carry_different_substituents_are_a_candidate():
    """N,N'-dimethyl-dimethylsulfodiimide: S(=NC)(=NC)(C)C, with one N-methyl grown to N-ethyl.

    Both pi neighbours are nitrogen, neutral, of the same isotope -- an atom-record comparison
    alone would call them interchangeable and drop the centre for good, because nothing later
    re-admits a candidate. So the local test may only reject a repeated pi direction when the
    neighbour is TERMINAL, where 'identical record' really does mean 'identical direction'.
    Everything else is the automorphism group's call.
    """
    m, sids = _mol(atoms='SNNCCCCC',
                   bonds=[(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1),
                          (1, 5, 1), (2, 6, 1), (6, 7, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['anchor'] == sids[0]


def test_alkyne_carbon_never_reaches_four_directions():
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (1, 2, 3), (2, 3, 1)])
    assert m.stereo_units() == []


def test_a_triple_bond_refuses_the_atom_even_when_the_count_would_reach_four():
    """A deliberately hypervalent phosphorus: P(#C)(C)(C)C, five bonds' worth of electrons.

    No real compound reaches four directions with a triple bond among them -- the valence is
    already spent -- so this is a forged record, in the spirit of test_derive.py's implicit_h=9.
    It exists because it is the only way to observe the early-out: count the triple bond as one
    direction and this record becomes a candidate.
    """
    m, sids = _mol(atoms='PCCCC', bonds=[(0, 1, 3), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
    assert m.stereo_units() == []


def test_an_over_coordinated_atom_does_not_overrun_the_direction_arrays():
    """Five, six, seven heavy neighbours, and sixteen drawn hydrogens. The point is the bounds
    check, not the chemistry.

    `refs` and the explicit-hydrogen array are four-element STACK arrays inside a `nogil`
    function, and the direction count that rejects a hypercoordinate record is only computed
    after the walk has finished writing into them. So the `n_ref < 4` and `n_h < 4` guards are
    the only thing between such a record and a stack buffer overflow.

    PF5 and SF6 are real compounds; the seven-fluorine phosphorus and the sixteen-hydrogen carbon
    are forged, in the spirit of
    `test_a_triple_bond_refuses_the_atom_even_when_the_count_would_reach_four`. Hydrogens are
    stated as zero throughout so that nothing here depends on an implicit count.

    How far past the array each case reaches is why those two counts are what they are. A short
    overrun only corrupts the frame's other locals, which the next atom overwrites anyway, so
    nothing observes it; a long one reaches the frame's own guard and the process dies. Measured
    with each guard removed in turn: seven heavy neighbours abort, and the hydrogen array needs
    twelve. Sixteen leaves margin, because where the compiler put the two arrays relative to each
    other is not something a test can pin.
    """
    pf5, _ = _mol(atoms='PFFFFF', hydrogens=[0] * 6,
                  bonds=[(0, k, 1) for k in range(1, 6)])
    assert pf5.stereo_units() == []

    sf6, _ = _mol(atoms='SFFFFFF', hydrogens=[0] * 7,
                  bonds=[(0, k, 1) for k in range(1, 7)])
    assert sf6.stereo_units() == []

    pf7, _ = _mol(atoms='PFFFFFFF', hydrogens=[0] * 8,
                  bonds=[(0, k, 1) for k in range(1, 8)])
    assert pf7.stereo_units() == []

    # the same bound applies to the hydrogen array, which ruling F26 gave its own four slots
    ch16, _ = _mol(atoms='C' + 'H' * 16, hydrogens=[0] * 17,
                   bonds=[(0, k, 1) for k in range(1, 17)])
    assert ch16.stereo_units() == []


def test_a_dative_bond_refuses_the_atom():
    """Trimethylamine-borane, N(C)(C)(C)->B, with the adduct bond as order 8.

    Three sigma carbons plus the dative bond would be four directions, and the nitrogen is
    four-coordinate exactly as an ammonium's is. Order 8 is refused anyway: it is chython's
    'anything else' order, carrying no geometry this rule could rely on. See the report --
    admitting it is a decision for the epic, not for this task.
    """
    m, sids = _mol(atoms='NCCCB', hydrogens=[0, 0, 0, 0, 3],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 8)])
    assert m.stereo_units() == []


def test_every_anchor_appears_at_most_once():
    """1,2-dibromo-1,2-difluoroethane: two centres, two units, two distinct anchors.

    Anchor uniqueness is what makes `unit_of` a function and what bounds the perception
    scratch at one record per atom; `test_the_anchor_collision_invariant_is_asserted` covers
    the enforcement, this covers the ordinary case.
    """
    m, sids = _mol(atoms='CFBrCFBr', hydrogens=[1, 0, 0, 1, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1)])
    units = m.stereo_units()
    anchors = [u['anchor'] for u in units]
    assert sorted(anchors) == sorted([sids[0], sids[3]])
    assert len(set(anchors)) == len(anchors)


def test_the_anchor_collision_invariant_is_asserted():
    """Two units on one anchor must raise, not be stored.

    Tetrahedral perception produces one kind and cannot collide, but the record has room for
    exactly one unit per anchor because a parity is keyed by anchor slot in SEG_PARITY. SU_CIS_TRANS
    (ruling F43) and SU_ATROPISOMER (ruling F45) anchor on atoms that already carry SU_TETRA in
    exactly the molecules those rulings were written to handle; the probe is what asserts the
    invariant their refusals depend on.
    """
    with pytest.raises(RuntimeError, match='anchor'):
        _core._stereo_anchor_collision_probe()


def test_a_second_read_is_idempotent_and_appends_nothing():
    """Two claims, and the second is the one that says the table is CACHED.

    Reading twice gives the same answer AND does not grow the arena. The first half alone would
    pass against an implementation that re-perceived every time -- it only rules out the second
    `structure_append` raising 'segment already attached' -- so the `total_len` half is what
    actually guards `ensure_stereo_units`' `structure_has` early return.

    The empty case is in here on purpose: an empty table is eight bytes and not zero, because
    `structure_append` of length 0 leaves the segment absent, `structure_has` keeps reporting it
    missing, and perception would re-run on every call.
    """
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    before = m.total_len
    first = m.stereo_units()
    grown = m.total_len
    assert grown > before
    assert m.stereo_units() == first
    assert m.total_len == grown

    empty, _ = _mol(atoms='CClCl', bonds=[(0, 1, 1), (0, 2, 1)])
    before = empty.total_len
    assert empty.stereo_units() == []
    grown = empty.total_len
    assert grown > before
    assert empty.stereo_units() == []
    assert empty.total_len == grown


def test_an_edit_rebuilds_the_table():
    """CHFCl has three directions; give it a bromine and it has four.

    The table is a cache over the arena, and an edit replaces the arena, so the cache cannot
    survive it. A stale one would report the wrong answer in both directions -- no unit on a
    centre that just became one, and a unit on one that stopped being one.
    """
    m, sids = _mol(atoms='CFCl', hydrogens=[1, 0, 0], bonds=[(0, 1, 1), (0, 2, 1)])
    assert m.stereo_units() == []
    with m.edit():
        br = m.add_atom('Br')
        m.add_bond(sids[0], br, 1)
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['refs'] == (sids[1], sids[2], br, None)
    with m.edit():
        m.delete_atom(br)
    assert m.stereo_units() == []


def test_remap_keeps_the_table_under_the_new_ids():
    """remap relabels ids without moving atom slots, and the table is keyed by slot.

    So the units survive the relabelling and come back under the new ids -- which is what makes
    a stored parity survive a remap too.
    """
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert len(m.stereo_units()) == 1
    m.remap({s: s + 100 for s in sids})
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['anchor'] == sids[0] + 100
    assert units[0]['refs'] == (sids[1] + 100, sids[2] + 100, sids[3] + 100, None)


def test_unit_of_finds_by_anchor():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.unit_of(sids[0])['kind'] == 0
    assert m.unit_of(sids[1]) is None


def test_unit_of_rejects_an_unknown_atom():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    with pytest.raises(KeyError):
        m.unit_of(max(sids) + 1)


def test_an_empty_molecule_has_no_units():
    m = MoleculeContainer()
    assert m.stereo_units() == []
    assert m.stereo_units() == []


def test_the_table_is_derived_and_does_not_travel_in_the_bytes():
    """to_bytes carries the persistent prefix only, and the unit table is not in it.

    Three claims. The packed length does not change when the table is present, and neither does
    the payload past the header. A structure rebuilt from those bytes re-derives the same table --
    that is what catches a stale segment entry surviving the round trip, because the header DOES
    travel inside the prefix and an entry left pointing past the end of the new buffer would be
    read as a live table.

    The comparison is over the WHOLE buffer, header included. It was written as `before[128:]` when
    v3's `total_len` and derived table entries lived inside the prefix and made a read change the
    header; v4 keeps derived segments out of the buffer, so the header is stable too and there is no
    reason to exempt it. Slicing at a hard-coded 128 would now also be wrong in a second way -- this
    molecule's header is 48 bytes, so the slice would have skipped three atom records.
    """
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    before = m.to_bytes()
    units = m.stereo_units()
    after = m.to_bytes()
    assert len(before) == len(after) == m.persistent_len < m.total_len
    assert before == after

    back = MoleculeContainer.from_bytes(after)
    assert back.stereo_units() == units


def test_the_segment_table_still_fits_the_locked_header():
    # 13 table entries * 8 bytes + 24 bytes of scalars == 128, the same size v3's header was, and
    # the table still starts at offset 24 -- so no persistent offset moved across the format change.
    # 128 is the C struct's size and the header's CEILING; this build writes 24 + 8 * seg_count and
    # stops the table one past its highest used segment, so most headers are 48 bytes.
    # SEG_STEREO_UNIT is absent from the table (it is derived, and derived caches have allocations of
    # their own), which is why its id may exceed the table width.
    #
    # THE TOTAL SEGMENT COUNT IS NOT ASSERTED HERE.  A literal count is a second copy of what
    # `test_structure.py` pins as the relation "the persistent ids are a dense prefix and the derived
    # ids follow", and appending one persistent segment renumbers the derived block -- so the literal
    # fails for a change it protects nothing against.  What this test is *for* is the header geometry.
    assert _core._header_size() == 128
    assert _core._segment_table_max() == 13
    assert _core.SEG_STEREO_UNIT >= _core._persistent_segment_count()
    # 4 bytes of scalars + anchor + 4 refs; the walks write this layout; `translate_stereo`,
    # `mark_stereogenic` and `_unit_dict` read it
    assert _core._stereo_unit_record_size() == 24


# ---------------------------------------------------------------------------------------------
# Cumulenes and atropisomers.
#
# Every ring below is written in KEKULE form, with alternating orders, because that is the only
# form the arena accepts: `add_bond(..., 4)` is journal-only and the apply raises
# NotImplementedError on it (`_molecule_container.pxi:428`), and `structure_from_bytes` rejects a
# stored HE_AROMATIC flag outright. So an aromatic ring here IS cyclohexa-1,3,5-triene, and the
# tests that say benzene has no cumulene unit are testing the ring rule, not an aromatic flag.


def test_but_2_ene_is_a_cis_trans_candidate():
    # CC=CC: a 2-atom chain, even, so cis/trans with the parity on the lower terminal
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[1] and u['n_refs'] == 4


def test_ethene_is_not_a_candidate():
    m, sids = _mol(atoms='CC', bonds=[(0, 1, 2)])
    assert m.stereo_units() == []


def test_penta_2_3_diene_is_an_axial_candidate():
    # CC=C=CC: a 3-atom chain, odd, so axial with the parity on the centre atom
    m, sids = _mol(atoms='CCCCC', bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 2 and u['anchor'] == sids[2]


def test_hexa_2_3_4_triene_is_a_cis_trans_candidate():
    # a 4-atom chain, even -- CT4 in the Blue Book's notation, parity on the central bond,
    # which this design stores on the lower terminal
    m, sids = _mol(atoms='CCCCCC',
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 2), (4, 5, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 1
    assert units[0]['anchor'] == sids[1]


def test_five_atom_cumulene_is_axial():
    m, sids = _mol(atoms='C' * 7,
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 2), (4, 5, 2), (5, 6, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 2 and units[0]['anchor'] == sids[3]


def test_symmetric_terminal_blocks_a_cumulene():
    # (CH3)2C=CHCH3 -- the left terminal's two directions are both methyl
    m, sids = _mol(atoms='CCCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 2), (3, 4, 1)])
    assert m.stereo_units() == []


def test_a_cumulene_terminal_orders_its_pair_by_ruling_f26_and_not_by_value():
    """A terminal's two directions: heavy slot first, then the hydrogen slot -- never sorted by
    raw numeric value.

    Both terminals here carry an explicit hydrogen declared BEFORE the heavy substituent, so the
    hydrogen's slot is the smaller number on the left terminal. Ruling F26 still puts the heavy
    direction first. Sorting the pair numerically -- which happens to look right whenever the
    absent direction is SU_NO_REF = 0xFFFFFFFF -- would emit (H, C) here and re-base the parity of
    every drawn-hydrogen cumulene relative to its implied-hydrogen twin.
    """
    m, sids = _mol(atoms='CHCCCH', hydrogens=[0] * 6,
                   bonds=[(0, 1, 1), (0, 2, 2), (0, 3, 1), (2, 4, 1), (2, 5, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[0]
    # lower terminal (slot 0) first: its heavy C then its H; then the other terminal's heavy C
    # then its H.  sids[3] before sids[1] is the whole point.
    assert u['refs'] == (sids[3], sids[1], sids[4], sids[5])
    assert u['unnamed_mask'] == 0


def test_an_implicit_hydrogen_on_a_cumulene_terminal_is_unnamed_and_sorts_last():
    """The same molecule with the hydrogens implied instead of drawn.

    The heavy directions keep their two positions and the hydrogen keeps its, so a parity stored
    against either record names the same configuration.  Only the names change.

    The mask is `0b1010` -- slots 1 and 3, one unnamed direction in each of the two pairs. Ruling
    F41's whole point is that this is not the same record as one with two unnamed directions in the
    same pair; `test_two_unnamed_directions_on_one_terminal_differ_from_one_on_each` is that pair.
    """
    m, sids = _mol(atoms='CCCC', hydrogens=[1, 0, 1, 0],
                   bonds=[(0, 1, 1), (0, 2, 2), (2, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['refs'] == (sids[1], None, sids[3], None)
    assert units[0]['unnamed_mask'] == 0b1010
    assert _popcount(units[0]['unnamed_mask']) == 2


def test_an_absent_direction_is_not_an_unnamed_one():
    """CC=CC forged with no hydrogens at all: each terminal has ONE direction, not two.

    `refs` is SU_NO_REF in both empty slots -- the entry means 'no atom here' whatever the reason
    -- but `unnamed_mask` marks implicit hydrogens, and there are none. Conflating the two would
    report both slots here and make the per-kind automorphism test (Ruling F34) unanswerable,
    since an empty slot is not a protium.
    """
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['refs'] == (sids[0], None, sids[3], None)
    assert units[0]['unnamed_mask'] == 0


def test_two_unnamed_directions_on_one_terminal_differ_from_one_on_each():
    """Propene against but-2-ene: Ruling F41's case, and why `spare` cannot be a count.

    Both records have two unnamed directions. In propene they are the =CH2 terminal's two implicit
    hydrogens, one direction twice over, and there is no cis/trans isomer; in but-2-ene there is one
    on each terminal and both configurations exist. A scalar `2` is the same answer for both, and
    `0b0011` against `0b1010` is not -- which is the whole of the ruling, because the automorphism
    filter decides stereogenicity from this field.

    Propene's pair never reaches a record, because `_terminal_pair` refuses a terminal whose every
    direction is unnamed (see `test_a_ch2_terminal_is_not_a_terminal`) -- so the two halves are
    asserted from the two sides available: propene has no unit at all, and but-2-ene's mask names the
    two SLOTS rather than counting them. A record with two unnamed directions in one pair is
    unreachable for exactly the reason propene is not stereogenic.
    """
    propene, _ = _mol(atoms='CCC', hydrogens=[3, 1, 2], bonds=[(0, 1, 1), (1, 2, 2)])
    assert all(u['kind'] != 1 for u in propene.stereo_units())

    # the methyl hydrogens are left off here, as elsewhere in this file, so that the methyl carbons
    # do not reach four directions and add tetrahedral records of their own
    butene, bsids = _mol(atoms='CCCC', hydrogens=[0, 1, 1, 0], bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    units = butene.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 1 and units[0]['anchor'] == bsids[1]
    assert units[0]['refs'] == (bsids[0], None, bsids[3], None)
    assert units[0]['unnamed_mask'] == 0b1010 and _popcount(units[0]['unnamed_mask']) == 2


def test_an_oxime_marks_the_hydrogen_slot_and_not_the_lone_pair_slot():
    """Acetaldoxime `CC=NO`: slot 1 is a real direction with no atom, slot 3 is no direction at all.

    This is the largest real population of E/Z units and the case a count cannot describe. The
    carbon terminal is (CH3, implicit H), so slot 1 is an implicit hydrogen; the nitrogen terminal is
    (O, -- ), and its second position is the lone pair, which `_terminal_pair` deliberately does not
    count as a direction. Both slots read SU_NO_REF in `refs`. Only the mask tells them apart, and
    that its bit 3 is CLEAR while bit 1 is set is Ruling F41 in one line.
    """
    m, sids = _mol(atoms='CCNO', hydrogens=[0, 1, 0, 1],
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 1 and u['anchor'] == sids[1]
    assert u['refs'] == (sids[0], None, sids[3], None)
    assert u['unnamed_mask'] == 0b0010 and _popcount(u['unnamed_mask']) == 1
    assert u['unnamed_mask'] & 0b1000 == 0


def test_a_sulfilimine_is_a_tetrahedral_centre_and_never_a_cumulene_terminal():
    """S-ethyl-S-methyl-N-methylsulfilimine, in both atom orders. Ruling F43.

    The sulfur has two sigma directions and no hydrogen, so it reads as a cumulene terminal, while
    the same atom reaches four directions in tetrahedral perception as two sigma, one pi and a
    lone pair and is already anchored SU_TETRA. Both units would key on it; a parity is keyed
    by anchor slot in SEG_PARITY, so one atom cannot hold two configurations.  Before the refusal
    this molecule raised out of `stereo_units()`, and only in the order where the sulfur holds
    the lower of the two terminal slots.  Both orders are here because only one of them is the
    regression.

    The refusal is on the terminal's own merits and not to dodge that collision: a cumulene
    terminal's two directions are its IN-PLANE SIGMA positions and a pyramidal sulfur has no plane
    for them to lie in. `R2S=NR` holds its configuration at the sulfur, which the tetrahedral record
    already names, so the putative E/Z across `S=N` would name nothing new.
    """
    # S=0: S(=N-CH3)(CH3)(CH2CH3).  The hydrogens are left off the carbons for the usual reason --
    # a stated CH3 or CH2 reaches four directions and adds a tetrahedral record of its own, which
    # would say nothing about this rule.
    sulfur_lower, a = _mol(atoms='SNCCCC', hydrogens=[0, 0, 0, 0, 0, 0],
                           bonds=[(0, 1, 2), (1, 2, 1), (0, 3, 1), (0, 4, 1), (4, 5, 1)])
    units = sulfur_lower.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 0 and units[0]['anchor'] == a[0]

    # the same molecule with the nitrogen declared first, which perceived cleanly even before F43
    # because the cis/trans anchor is the LOWER terminal and that was then the nitrogen
    nitrogen_lower, b = _mol(atoms='NSCCCC', hydrogens=[0, 0, 0, 0, 0, 0],
                             bonds=[(0, 1, 2), (0, 2, 1), (1, 3, 1), (1, 4, 1), (4, 5, 1)])
    units = nitrogen_lower.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 0 and units[0]['anchor'] == b[1]


def test_a_cumulene_terminal_with_three_substituents_is_refused():
    """A forged sp2 terminal with three sigma neighbours besides the chain: not a terminal.

    Four in-plane directions is not a geometry this record can describe, and admitting it would
    silently keep only the first two. The atom is still a tetrahedral candidate -- three sigma
    neighbours plus one pi is the tetrahedral perception's four directions -- which is why this
    asserts on the kind rather than on an empty list.
    """
    m, sids = _mol(atoms='CFFFCC', hydrogens=[0] * 6,
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 2), (4, 5, 1)])
    units = m.stereo_units()
    assert all(u['kind'] != 1 and u['kind'] != 2 for u in units)
    assert [u['anchor'] for u in units] == [sids[0]]


def test_a_triple_bond_breaks_a_cumulene_chain():
    # CC=C#CC: the chain stops at the triple bond, and a terminal carrying one is refused
    m, sids = _mol(atoms='CCCCC', bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 4, 1)])
    assert m.stereo_units() == []


# Kekule ring bond lists, named the same way `_BIPHENYL_BONDS` below is, because each of them is
# used by more than one test and every one of them has to alternate correctly by hand -- an order-4
# bond raises at apply, so there is no aromatic spelling available here.
_BENZENE_BONDS = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)]
# naphthalene: atoms 0 and 9 are the fusion pair, 1..8 the peripheral carbons
_NAPHTHALENE_BONDS = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 9, 2),
                      (9, 8, 1), (8, 7, 2), (7, 6, 1), (6, 5, 2), (5, 0, 1),
                      (9, 0, 1)]


def test_benzene_has_no_cumulene_unit():
    """Kekule benzene is cyclohexa-1,3,5-triene in the arena, and none of its three double bonds
    is a cis/trans candidate.

    This is Ruling F31's observable half. The trap it names -- a chain walk keyed on `order == 2`
    running through an aromatic ring -- takes a different shape than the ruling assumed, because
    kekulisation ALTERNATES: benzene's order-2 half-edges form three separate two-atom chains, not
    one six-atom chain. Each of those three chains still passes the terminal test (one ring
    neighbour and one hydrogen, plainly distinguishable), so without the ring rule benzene reports
    three cis/trans units and toluene four. `HE_AROMATIC` cannot be the discriminator on this
    branch: nothing sets it, `structure_from_bytes` rejects it, and an order-4 bond raises on
    apply, so a Kekule aromatic ring and a real cyclohexatriene are the same graph here.
    """
    m, sids = _mol(atoms='C' * 6, hydrogens=[1] * 6, bonds=_BENZENE_BONDS)
    assert m.stereo_units() == []


def test_naphthalene_has_no_cumulene_unit():
    m, sids = _mol(atoms='C' * 10, hydrogens=[0, 0, 1, 1, 1, 1, 1, 1, 0, 0],
                   bonds=_NAPHTHALENE_BONDS)
    assert all(u['kind'] != 1 and u['kind'] != 2 for u in m.stereo_units())


def test_toluene_has_no_cumulene_unit_only_its_methyl_carbon():
    """The same molecule `test_the_unnamed_direction_slots_are_recorded` pins, asserted from
    the other side.

    `test_the_unnamed_direction_slots_are_recorded` says toluene has exactly one unit; drop
    the ring rule and it has four, so that test fails too. This one says WHICH kinds are absent, so
    the failure names the cause.
    """
    m, sids = _mol(atoms='CCCCCCC', hydrogens=[3, 0, 1, 1, 1, 1, 1],
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2),
                          (4, 5, 1), (5, 6, 2), (6, 1, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 0 and units[0]['anchor'] == sids[0]


def test_cyclohexene_and_cyclooctene_split_on_the_ring_threshold():
    """The threshold is a ring smaller than 8, so this pair is what makes it observable.

    Cyclohexene's double bond cannot be trans -- the ring holds its two substituents cis -- while
    trans-cyclooctene is a real, isolable compound. Excluding every ring double bond instead of
    small-ring ones would lose the second, and with it the epic's macrocyclic-cumulene gate.
    """
    hexene = [(0, 1, 2)] + [(k, k + 1, 1) for k in range(1, 5)] + [(5, 0, 1)]
    m6, s6 = _mol(atoms='C' * 6, hydrogens=[1, 1, 2, 2, 2, 2], bonds=hexene)
    assert all(u['kind'] != 1 for u in m6.stereo_units())

    octene = [(0, 1, 2)] + [(k, k + 1, 1) for k in range(1, 7)] + [(7, 0, 1)]
    m8, s8 = _mol(atoms='C' * 8, hydrogens=[1, 1] + [2] * 6, bonds=octene)
    units = [u for u in m8.stereo_units() if u['kind'] == 1]
    assert len(units) == 1 and units[0]['anchor'] == s8[0]


def test_a_macrocyclic_allene_survives_the_ring_rule():
    """A 1,2-cyclotridecadiene: an axial cumulene inside a 13-ring.

    Ring cumulenes in rings of 8 and up are exactly what the epic's gate 3 pins, so the ring rule
    must not reach them.
    """
    bonds = [(0, 1, 2), (1, 2, 2)] + [(k, k + 1, 1) for k in range(2, 12)] + [(12, 0, 1)]
    m, sids = _mol(atoms='C' * 13, hydrogens=[1, 0, 1] + [2] * 10, bonds=bonds)
    units = [u for u in m.stereo_units() if u['kind'] == 2]
    assert len(units) == 1 and units[0]['anchor'] == sids[1]


def test_a_small_ring_allene_keeps_its_axial_candidate():
    """1,2-cyclohexadiene and 1,2-cycloheptadiene: Ruling F44.

    The small-ring cut's argument -- the ring path holds the terminals' substituents cis, so there is
    no second configuration for a parity to name -- is a statement about a CIS/TRANS unit and is
    meaningless for an axial one. An allene's terminals are perpendicular, 'cis' is not defined for
    them, and its two configurations are enantiomers that no ring path can equate. Applying the cut
    to both kinds lost the axial candidates of these two while
    `test_a_macrocyclic_allene_survives_the_ring_rule` kept the thirteen-ring's -- wrong on the
    strained members of the exact axis this epic is staked on. Both are isolable only in trapping
    experiments and both are chiral, with a literature on enantioselective capture.
    """
    six = [(0, 1, 2), (1, 2, 2)] + [(k, k + 1, 1) for k in range(2, 5)] + [(5, 0, 1)]
    m6, s6 = _mol(atoms='C' * 6, hydrogens=[1, 0, 1, 2, 2, 2], bonds=six)
    units = [u for u in m6.stereo_units() if u['kind'] == 2]
    assert len(units) == 1 and units[0]['anchor'] == s6[1]

    seven = [(0, 1, 2), (1, 2, 2)] + [(k, k + 1, 1) for k in range(2, 6)] + [(6, 0, 1)]
    m7, s7 = _mol(atoms='C' * 7, hydrogens=[1, 0, 1, 2, 2, 2, 2], bonds=seven)
    units = [u for u in m7.stereo_units() if u['kind'] == 2]
    assert len(units) == 1 and units[0]['anchor'] == s7[1]


def test_a_cyclopropane_fused_at_each_terminal_does_not_cut_a_twelve_ring_double_bond():
    """A twelve-ring double bond with a cyclopropane fused at each of its two terminals.

    Both terminals are on a three-ring and the two of them share the twelve-ring, so asking 'is
    either atom on a small ring' and 'do they share some ring' as INDEPENDENT questions refuses this
    bond -- while the cyclopropanes constrain nothing whatever about the twelve-ring, whose double
    bond has both configurations. The test has to be one test: intersect the two atoms' prototypes
    first, then measure the size of a prototype they actually SHARE.

    Losing a candidate is the expensive direction here, since nothing after perception re-examines a
    refused one.
    """
    bonds = [(0, 1, 2)] + [(k, k + 1, 1) for k in range(1, 11)] + [(11, 0, 1)]
    bonds += [(0, 12, 1), (12, 11, 1),      # three-ring fused on the 0-11 edge
              (1, 13, 1), (13, 2, 1)]       # three-ring fused on the 1-2 edge
    m, sids = _mol(atoms='C' * 14, hydrogens=[0, 0] + [2] * 12, bonds=bonds)
    units = [u for u in m.stereo_units() if u['kind'] == 1]
    assert [u['anchor'] for u in units] == [sids[0]]


_BIPHENYL_BONDS = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                   (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1),
                   (0, 6, 1)]
_BIPHENYL_H = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1]


def test_biphenyl_with_ortho_substituents_is_an_atropisomer():
    # 2,2'-dichloro/fluoro-biphenyl: two Kekule rings, a single acyclic pivot bond, ortho Cl and F
    m, sids = _mol(atoms='C' * 12 + 'ClF', hydrogens=_BIPHENYL_H + [0, 0],
                   bonds=_BIPHENYL_BONDS + [(1, 12, 1), (7, 13, 1)])
    units = m.stereo_units()
    assert len(units) == 1
    u = units[0]
    assert u['kind'] == 3 and u['anchor'] == sids[0] and u['n_refs'] == 4
    # the two ortho ring directions on each end, lower-slot pivot first, each pair CSR ascending
    assert u['refs'] == (sids[1], sids[5], sids[7], sids[11])
    assert u['unnamed_mask'] == 0


def test_unsubstituted_biphenyl_is_not_an_atropisomer():
    m, sids = _mol(atoms='C' * 12, hydrogens=[0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
                   bonds=_BIPHENYL_BONDS)
    assert m.stereo_units() == []


def test_an_explicit_ortho_hydrogen_is_not_an_ortho_substituent():
    """The same biphenyl with its two ortho hydrogens DRAWN instead of implied.

    A hydrogen is not what hinders rotation, and explicitness is a drawing choice, so counting a
    drawn one as the ortho substituent would make this molecule an atropisomer and its
    implicit-hydrogen twin above not one -- the same representation dependence Ruling F26 exists to
    keep out of the reference order.
    """
    m, sids = _mol(atoms='C' * 12 + 'HH', hydrogens=_BIPHENYL_H + [0, 0],
                   bonds=_BIPHENYL_BONDS + [(1, 12, 1), (7, 13, 1)])
    assert m.stereo_units() == []


def test_one_ortho_substituent_on_each_end_is_required():
    """Only one ring substituted: 2-chlorobiphenyl rotates freely enough to be one compound."""
    m, sids = _mol(atoms='C' * 12 + 'Cl', hydrogens=_BIPHENYL_H[:6] + [0, 1] + _BIPHENYL_H[8:] + [0],
                   bonds=_BIPHENYL_BONDS + [(1, 12, 1)])
    assert all(u['kind'] != 3 for u in m.stereo_units())


def test_ring_bond_is_never_an_atropisomer_axis():
    # decalin's fusion bond is in a ring, so it is excluded by rule
    m, sids = _mol(atoms='C' * 10,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                          (0, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 5, 1)])
    assert all(u['kind'] != 3 for u in m.stereo_units())


def test_a_fused_ring_bond_with_peri_substituents_is_not_an_axis():
    """1,8-dichloronaphthalene: the fusion bond passes every other clause of the rule.

    Its two atoms are ring atoms of degree three with no hydrogen, all their remaining bonds are
    ring bonds, and each end has a substituted neighbour -- so `HE_IN_RING` on the pivot bond is the
    ONLY clause that rejects it. Decalin does not isolate that clause, because its fusion bond is
    also rejected for having no ortho substituent anywhere. There is no axis here to be configured:
    the two rings are one rigid plane.
    """
    m, sids = _mol(atoms='C' * 10 + 'ClCl',
                   hydrogens=[0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                   bonds=_NAPHTHALENE_BONDS + [(1, 10, 1), (8, 11, 1)])
    assert all(u['kind'] != 3 for u in m.stereo_units())


def test_a_saturated_pivot_is_a_tetrahedral_centre_and_not_an_axis():
    """2,2'-dimethyl-bicyclohexyl: the collision case the anchor invariant would die on.

    A cyclohexane-cyclohexane pivot carbon has two ring bonds, the pivot bond and an implicit
    hydrogen -- four directions, so TETRAHEDRAL is already anchored there. If the
    atropisomer rule also fired on this axis it would claim the same anchor; a parity is keyed by
    anchor slot in SEG_PARITY, so one atom cannot hold two configurations, and perception would raise
    on an ordinary molecule.
    Requiring the pivot to carry no hydrogen direction is what keeps them apart, and it is also the
    right chemistry: a saturated C-C bond rotates, however crowded its ortho positions are.
    """
    bonds = ([(k, (k + 1) % 6, 1) for k in range(6)]
             + [(6 + k, 6 + (k + 1) % 6, 1) for k in range(6)]
             + [(0, 6, 1), (1, 12, 1), (7, 13, 1)])
    m, sids = _mol(atoms='C' * 14,
                   hydrogens=[1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3], bonds=bonds)
    units = m.stereo_units()
    assert all(u['kind'] != 3 for u in units)
    assert m.unit_of(sids[0])['kind'] == 0
    assert m.unit_of(sids[6])['kind'] == 0


def test_a_symmetric_ortho_pair_is_still_a_candidate_here():
    """2,6-dichloro-2'-chlorobiphenyl. Pinning what perception does NOT decide, not endorsing it.

    Both ortho positions of the left ring carry a chlorine, so rotating that ring by half a turn
    reproduces the molecule and the axis has no second configuration -- this is not an atropisomer.
    The rule's distinguishability clause does not catch it: that clause uses the same comparator the
    automorphism group uses, and the comparator answers only for TERMINAL atoms, while a pivot's
    ring neighbours have degree two or more. Deciding the two ortho directions apart needs the rings
    walked, which is stereogenicity -- the automorphism filter's job.  So this asserts that
    perception emits the candidate without deciding it; `mark_stereogenic` decides stereogenicity,
    and for this molecule it finds the unit stereogenic (the right ring's asymmetry breaks the
    symmetry of the left ring's ortho pair).
    """
    hydrogens = list(_BIPHENYL_H)
    hydrogens[5] = 0
    m, sids = _mol(atoms='C' * 12 + 'ClClCl', hydrogens=hydrogens + [0, 0, 0],
                   bonds=_BIPHENYL_BONDS + [(1, 12, 1), (5, 13, 1), (7, 14, 1)])
    units = m.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 3 and units[0]['anchor'] == sids[0]


# Kekule cyclooctatetraene, twice: an eight-ring is the smallest aryl-like ring the small-ring cut
# does NOT reach, so a biaryl built from two of them is where the pivot's own ring double bond and
# the atropisomer axis want the same anchor.  Ring 1 is atoms 0-7, ring 2 atoms 8-15, and both lists
# alternate from the lower-numbered atom of each ring.
_COT_ONE = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1), (6, 7, 2), (7, 0, 1)]
_COT_TWO = [(8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 12, 1), (12, 13, 2), (13, 14, 1),
            (14, 15, 2), (15, 8, 1)]


def test_a_biaryl_of_eight_rings_relocates_its_atropisomer_anchor():
    """An ortho-dichloro bi(cyclooctatetraenyl), in both atom orders. Ruling F45.

    A Kekule aryl pivot ALWAYS carries a ring double bond, so it is always a cis/trans terminal as
    well as an axis end. Plain biphenyl never collides only because the small-ring cut drops a
    six-ring's chain; at ring size 8 the cut does not reach it and both units qualify, so a fixed
    anchor choice raises out of `stereo_units()` here. Refusing a pivot that has a ring double bond
    kills every biaryl instead.

    A bond kind may anchor at EITHER end, so the answer is anchor choice: the axis takes whichever
    pivot is free. Both orders are here because which pivot the cumulene pass claims depends on slot order
    -- in the first record it is the lower pivot, so the axis relocates; in the second it is the upper
    one, so the axis stays where it started. The record's MEANING does not move either way, because
    the parity is stored against `refs` whose two pairs follow the anchor.
    """
    # ring 1's pivot is atom 0, whose ring double bond makes it the lower terminal of that chain, so
    # the axis has to move to atom 15 -- ring 2's pivot, whose own double bond runs to atom 14
    lower_taken, a = _mol(atoms='C' * 16 + 'ClCl',
                          hydrogens=[0, 0] + [1] * 6 + [1] * 6 + [0, 0] + [0, 0],
                          bonds=_COT_ONE + _COT_TWO + [(0, 15, 1), (1, 16, 1), (14, 17, 1)])
    units = lower_taken.stereo_units()
    anchors = [u['anchor'] for u in units]
    assert len(anchors) == len(set(anchors)) == 9
    assert lower_taken.unit_of(a[0])['kind'] == 1        # the pivot's own ring double bond
    atrop = [u for u in units if u['kind'] == 3]
    assert len(atrop) == 1 and atrop[0]['anchor'] == a[15]
    assert atrop[0]['refs'] == (a[8], a[14], a[1], a[7])  # the anchor's ortho pair leads

    # the same shape with each ring turned so that the axis' UPPER pivot is the taken one: no
    # relocation happens and the axis anchors on the lower pivot, as every bond kind does by default
    upper_taken, b = _mol(atoms='C' * 16 + 'ClCl',
                          hydrogens=[1] * 6 + [0, 0, 0, 0] + [1] * 6 + [0, 0],
                          bonds=_COT_ONE + _COT_TWO + [(7, 8, 1), (6, 16, 1), (9, 17, 1)])
    units = upper_taken.stereo_units()
    anchors = [u['anchor'] for u in units]
    assert len(anchors) == len(set(anchors)) == 9
    assert upper_taken.unit_of(b[8])['kind'] == 1
    atrop = [u for u in units if u['kind'] == 3]
    assert len(atrop) == 1 and atrop[0]['anchor'] == b[7]

    # and the ordinary case is untouched: biphenyl's axis still anchors on its lower pivot, which is
    # free because the six-ring's chain is cut
    biphenyl, c = _mol(atoms='C' * 12 + 'ClF', hydrogens=_BIPHENYL_H + [0, 0],
                       bonds=_BIPHENYL_BONDS + [(1, 12, 1), (7, 13, 1)])
    units = biphenyl.stereo_units()
    assert len(units) == 1 and units[0]['kind'] == 3 and units[0]['anchor'] == c[0]


def test_a_biaryl_of_eight_rings_can_lose_its_axis_to_both_pivots():
    """The other side of ruling F45's relocation: when BOTH pivots are taken, the axis is dropped.

    The relocation above needs one free pivot. Join the two eight-rings at the atom each ring's
    Kekule chain starts from -- atoms 0 and 8 -- and each pivot is the LOWER terminal of its own ring
    double bond, so the cumulene pass anchors at both and the axis has nowhere to go. The shape is
    REACHABLE, which is why this test exists. Both bond orientations of the rings are here because the
    answer must not depend on which way round the Kekule chain is written.

    What is asserted is the CURRENT behaviour, not the ideal one: perception is allowed to lose a
    candidate, and the principled fix (relocate the colliding cis/trans unit to its own other
    terminal) is a cascade that this branch does not build. If a later task builds it, this test
    changes -- deliberately, with the axis appearing on one of the two pivots.
    """
    hydrogens = [0, 0] + [1] * 6 + [0, 0] + [1] * 6 + [0, 0]
    for chain_from_the_pivot in (True, False):
        one, two = (_COT_ONE, _COT_TWO) if chain_from_the_pivot else (
            [(i, j, 3 - o) for i, j, o in _COT_ONE], [(i, j, 3 - o) for i, j, o in _COT_TWO])
        m, a = _mol(atoms='C' * 16 + 'ClCl', hydrogens=hydrogens,
                    bonds=one + two + [(0, 8, 1), (1, 16, 1), (9, 17, 1)])
        units = m.stereo_units()
        # every ring double bond is a cis/trans unit and nothing else is emitted
        assert {u['kind'] for u in units} == {1}
        assert len(units) == 8
        # both pivots are among the anchors -- that is precisely why the axis was refused
        assert a[0] in {u['anchor'] for u in units} and a[8] in {u['anchor'] for u in units}


# 1,1'-binaphthyl: two `_NAPHTHALENE_BONDS` rings, the second offset by ten, joined at the atoms
# numbered 1 and 11 -- each of which is a C1 position, ortho to the peri fusion atom of its own ring.
_BINAPHTHYL_BONDS = (_NAPHTHALENE_BONDS
                     + [(i + 10, j + 10, o) for i, j, o in _NAPHTHALENE_BONDS]
                     + [(1, 11, 1)])
# atoms 0 and 9 are the fusion pair and atom 1 is the pivot, so those three carry no hydrogen; the
# rest do, and atom 5's is the PERI hydrogen that is the whole of unsubstituted binaphthyl's barrier
_BINAPHTHYL_H = [0, 0, 1, 1, 1, 1, 1, 1, 1, 0] * 2


def test_binaphthyl_is_an_atropisomer_through_its_peri_fusion_atom():
    """1,1'-binaphthyl and BINOL. Ruling F46.

    Unsubstituted binaphthyl's rotational barrier is the PERI hydrogen on C8, so nothing hangs off
    the fusion atom C8a as an exocyclic substituent and the plain ortho test -- which looks only at
    bonds leaving the ring -- perceived nothing at all, while BINOL, the same scaffold plus two
    hydroxyls, was already a candidate. Binaphthyl and BINAP are most of why this kind exists, so
    half of it was invisible. A ring-fusion ortho neighbour therefore counts as hindering.

    The two are asserted together because they are the discrimination: remove the fusion rule and
    BINOL keeps its unit and binaphthyl loses its.
    """
    m, sids = _mol(atoms='C' * 20, hydrogens=_BINAPHTHYL_H, bonds=_BINAPHTHYL_BONDS)
    units = m.stereo_units()
    assert len(units) == 1
    assert units[0]['kind'] == 3 and units[0]['anchor'] == sids[1]
    # the anchor's ortho pair leads: the fusion atom and the CH, CSR ascending
    assert units[0]['refs'] == (sids[0], sids[2], sids[10], sids[12])

    # BINOL: 1,1'-bi-2-naphthol, hydroxyls on the two carbons ortho to the axis on the other side.
    # Those two carbons lose the hydrogen they had, or they reach four directions and add
    # tetrahedral records that say nothing about this rule.
    hydrogens = list(_BINAPHTHYL_H)
    hydrogens[2] = hydrogens[12] = 0
    binol, bsids = _mol(atoms='C' * 20 + 'OO', hydrogens=hydrogens + [1, 1],
                        bonds=_BINAPHTHYL_BONDS + [(2, 20, 1), (12, 21, 1)])
    units = binol.stereo_units()
    assert len(units) == 1
    assert units[0]['kind'] == 3 and units[0]['anchor'] == bsids[1]


def test_every_anchor_is_distinct():
    """A molecule with one of each kind: no atom may anchor two units.

    CHFClBr, but-2-ene and penta-2,3-diene in one record -- a tetrahedral centre, an even cumulene
    and an odd one. Three units, three anchors, and `unit_of` is a function of the atom only
    because that holds.
    """
    m, sids = _mol(atoms='CFClBr' + 'CCCC' + 'CCCCC', hydrogens=[1] + [0] * 12,
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1),
                          (4, 5, 1), (5, 6, 2), (6, 7, 1),
                          (8, 9, 1), (9, 10, 2), (10, 11, 2), (11, 12, 1)])
    units = m.stereo_units()
    anchors = [u['anchor'] for u in units]
    assert len(anchors) == len(set(anchors)) == 3
    assert sorted(u['kind'] for u in units) == [0, 1, 2]
    assert m.unit_of(sids[0])['kind'] == 0
    assert m.unit_of(sids[5])['kind'] == 1
    assert m.unit_of(sids[10])['kind'] == 2


def test_the_single_gate_refuses_a_second_unit_on_one_anchor():
    """`_stereo_emit` is the ONLY enforcement of the anchor no-collision invariant (Ruling F42).

    A second gate over the finished table cannot fire while `_stereo_emit` is the only emitter, so the
    gate is here and the probe exercises the gate itself rather than a helper beside it.

    No molecule reaches it, which is the point -- every kind that could collide either refuses
    locally (a sulfur cumulene terminal, Ruling F43) or relocates its anchor (an atropisomer whose
    lower pivot is taken, Ruling F45). A probe is therefore the only way the live gate is observable.
    """
    with pytest.raises(RuntimeError, match='anchor'):
        _core._stereo_anchor_collision_probe()


def test_a_cumulene_terminal_does_not_overrun_its_two_ref_slots():
    """A forged terminal with twenty-four sigma neighbours, and one with twenty-four drawn
    hydrogens.

    Each terminal of a cumulene writes its two directions into a two-word window of the unit's
    four-word `refs`, and the count that rejects an over-substituted terminal is only computed
    after the walk has finished writing. So the two `< 2` guards are the only thing between this
    record and a stack buffer overflow, and a short overrun is not observable -- Ruling F32:
    `test_an_over_coordinated_atom_does_not_overrun_the_direction_arrays` measured twelve drawn
    hydrogens before the hydrogen-array overrun reached the frame guard. Twenty-four leaves
    margin, since where the compiler puts the record relative to the frame's guard is not
    something a test can pin.

    Both records are forged, in the spirit of
    `test_an_over_coordinated_atom_does_not_overrun_the_direction_arrays`.
    """
    heavy, _ = _mol(atoms='C' + 'F' * 24 + 'C', hydrogens=[0] * 26,
                    bonds=[(0, k, 1) for k in range(1, 25)] + [(0, 25, 2)])
    assert heavy.stereo_units() == []

    drawn, _ = _mol(atoms='C' + 'H' * 24 + 'C', hydrogens=[0] * 26,
                    bonds=[(0, k, 1) for k in range(1, 25)] + [(0, 25, 2)])
    assert drawn.stereo_units() == []


def test_a_fully_cumulated_ring_terminates():
    """A forged four-ring of nothing but double bonds: every atom is a chain INTERIOR.

    The walk finds a chain by starting from an atom with exactly one double bond, so a cycle in the
    double-bond subgraph has no entry point and is never walked. An implementation that started
    anywhere and followed double bonds would loop forever here.
    """
    m, sids = _mol(atoms='C' * 4, hydrogens=[0] * 4,
                   bonds=[(0, 1, 2), (1, 2, 2), (2, 3, 2), (3, 0, 2)])
    assert m.stereo_units() == []


def test_a_ch2_terminal_is_not_a_terminal():
    """Propene and isobutene: the =CH2 end's two directions are both implicit protium.

    Neither has a cis/trans isomer, and the reason is the same as for two sulfur lone pairs --
    two implicit hydrogens are one direction twice over, and unlike a pair of heavy neighbours
    they can never become distinguishable later, so refusing here loses nothing.
    """
    propene, _ = _mol(atoms='CCC', hydrogens=[3, 1, 2], bonds=[(0, 1, 1), (1, 2, 2)])
    assert all(u['kind'] != 1 for u in propene.stereo_units())

    isobutene, _ = _mol(atoms='CCCC', hydrogens=[3, 0, 3, 2],
                        bonds=[(0, 1, 1), (1, 2, 1), (1, 3, 2)])
    assert all(u['kind'] != 1 for u in isobutene.stereo_units())


def test_a_branched_double_bond_subgraph_is_refused_and_does_not_hang():
    """A forged triangle of double bonds with a fourth double bond hanging off it.

    Atom 3 is the only atom with one chain bond, so the walk starts there and steps onto atom 0,
    whose chain degree is three. Without the refusal the walk goes 3, 0, 1, 2, 0, 1, 2, ... forever:
    at atom 0 the lowest-numbered chain neighbour that is not where it came from is always inside
    the triangle, so it never finds its way back out. Refusing a chain vertex of degree above two is
    the whole termination argument -- there is no step counter -- and this record is what tests it.
    The bond orders here are impossible for carbon; it is forged, in the spirit of
    `test_a_triple_bond_refuses_the_atom_even_when_the_count_would_reach_four`.
    """
    m, sids = _mol(atoms='C' * 4, hydrogens=[0] * 4,
                   bonds=[(0, 1, 2), (1, 2, 2), (2, 0, 2), (0, 3, 2)])
    assert all(u['kind'] != 1 and u['kind'] != 2 for u in m.stereo_units())


def test_a_cumulene_interior_with_a_third_neighbour_is_refused():
    """A forged allene centre carrying a methyl: an sp carbon has exactly two neighbours.

    Both terminals are deliberately well formed -- a methyl and an implicit hydrogen each, so each
    would pass `_terminal_pair` -- because a terminal that fails on its own would mask the interior
    rule and leave it untested. The only thing wrong with this record is the third bond on atom 2,
    and a real 2-methyl-penta-2,3-diene cannot exist: sp carbon has no third direction.
    """
    m, sids = _mol(atoms='CCCCCC', hydrogens=[3, 1, 0, 1, 3, 3],
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 1), (2, 5, 1)])
    assert all(u['kind'] != 2 for u in m.stereo_units())
