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
"""Implicit hydrogens: the authority order, and the ``vvv`` measurements behind it (written up in
``_hydrogens.py``).  A CTfile's stated total valence does *not* outrank the valence rules.  A count
nothing determines is unknown rather than zero and is never grounds for rejecting the record; the state
is the arena's ``H_UNKNOWN``, so ``implicit_h_of`` answers ``None`` and the atom is in
``unknown_hydrogens``.
"""

from pytest import skip

from ....core._core import (HYD_AMBIGUOUS_AROMATIC, HYD_DERIVED, HYD_DERIVED_OTHER_READING,
                            HYD_NO_AROMATIC_FORM, HYD_NO_VALENCE_RULE, HYD_REASON_MASK)
from .._hydrogens import (H_MAX, MRV_IMPLICIT_H, calc_implicit, implicit_for_atom,
                          implicit_from_valence, implicit_h_records, valence_for_write)
from .._sdf import parse_record
from .._sgroup import UNSUPPORTED
from .._v2000 import emit_v2000, parse_v2000
from .conftest import holds_an_aromatic_bond


def _record(atoms, bonds, extra=(), title='t'):
    """A V2000 record from ``[(symbol, vvv), ...]`` and ``[(a, b, order), ...]``, 1-based bonds."""
    lines = [title, '  test', '', f'{len(atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000']
    for i, (symbol, valence) in enumerate(atoms):
        v = 0 if valence is None else valence
        lines.append(f'{float(i):10.4f}    0.0000    0.0000 {symbol:<3} 0  0  0  0  0{v:3d}'
                     f'  0  0  0  0  0')
    for a, b, order in bonds:
        lines.append(f'{a:3d}{b:3d}{order:3d}  0  0  0  0')
    lines.extend(extra)
    lines.append('M  END')
    return lines


def _read(atoms, bonds, extra=()):
    ctab = parse_v2000(_record(atoms, bonds, extra), [])
    mol, store, log = ctab.build()
    return ctab, mol, log


def _counts(mol):
    return [mol.implicit_h_of(s) for s in mol.atom_numbers]


# ------------------------------------------------------------------------------- the ordinary case

def test_the_valence_rules_answer_when_nothing_is_stated():
    _, mol, log = _read([('C', None), ('O', None)], [(1, 2, 1)])
    assert _counts(mol) == [3, 1]
    assert not log


def test_a_charge_changes_the_answer():
    """The rule key includes the charge; an ammonium nitrogen carries four hydrogens, not three."""
    _, mol, _ = _read([('N', None)], [], ['M  CHG  1   1   1'])
    assert _counts(mol) == [4]


def test_an_order_eight_bond_contributes_nothing_to_the_valence_sum():
    """chython's order 8 is a dative or ionic contact and carries no electron pair.  Counting it
    would make the iron of every ferrocene a valence error and take a hydrogen off its partner."""
    _, mol, _ = _read([('C', None), ('N', None)], [(1, 2, 8)])
    assert _counts(mol)[0] == 4, 'the carbon still has four hydrogens; the contact is not a bond order'


# --------------------------------------------------------- the authority order, and why it is that

def test_a_stated_hydrogen_count_outranks_everything():
    """``MRV_IMPLICIT_H`` is the one channel in either MDL version that states a count as a count."""
    extra = ['M  STY  1   1 DAT', 'M  SAL   1  1   1', f'M  SDT   1 {MRV_IMPLICIT_H}',
             'M  SED   1 IMPL_H1']
    _, mol, _ = _read([('C', None)], [], extra)
    assert _counts(mol) == [1], 'the file said one hydrogen; nothing here second-guesses it'


def test_the_stated_count_is_found_whatever_case_its_field_name_is_written_in():
    """The field name is an extension's keyword, so another spelling of it is the same statement; matching
    it exactly falls back to the valence rules on the one atom the field exists to rescue.  The fold is
    *this field's* only -- a ``FIELDNAME`` is free text in the spec, where case can be the whole meaning.
    """
    extra = ['M  STY  1   1 DAT', 'M  SAL   1  1   1', 'M  SDT   1 mrv_implicit_h',
             'M  SED   1 IMPL_H1']
    _, mol, log = _read([('C', None)], [], extra)
    assert _counts(mol) == [1]
    assert any('same field' in x for x in log), log


def test_a_stated_valence_does_not_override_the_valence_rules():
    """The ``vvv`` measurement: of the 42 corpus atoms stating a valence, ``vvv - drawn`` disagrees with the
    rules on 21 and is arithmetically impossible on 2 more, and the rules are right wherever chemistry can
    check.  Two ferrocene records in ``test/standardize.sdf`` show why -- record 76 draws the Kekule double
    bonds, record 78 carries the same ``vvv`` on an all-single ring -- so ``vvv`` counts orders the file
    does not always draw, and the subtraction is a hydrogen count only when the drawing is complete.
    """
    # a Cp carbon as record 78 draws it: two single bonds, vvv=4.  Honouring vvv gives 2 hydrogens.
    _, mol, log = _read([('C', 4), ('C', None), ('C', None)], [(1, 2, 1), (1, 3, 1)])
    assert _counts(mol)[0] == 2  # sic: with two single bonds drawn the rules also say 2
    # and now the case that separates them -- a fully drawn aromatic carbon.
    _, mol, log = _read([('C', 4), ('C', None), ('C', None)], [(1, 2, 2), (1, 3, 1)])
    assert _counts(mol)[0] == 1, 'the rules give 1; vvv - drawn would give 1 here too'
    # the real divergence: a valence stating more than the drawing accounts for.
    _, mol, log = _read([('C', 6), ('C', None)], [(1, 2, 1)])
    assert _counts(mol)[0] == 3, 'the rules give 3; honouring vvv=6 would give 5'
    assert any('stated valence 6' in x and 'the derivation gives 3' in x and 'using 3' in x
               for x in log), log
    # the line names the derivation and not "the valence rules": an aromatic atom's count comes from the
    # kekuliser's classifier and then a row, so naming the row alone sends a reader to the wrong table.
    assert not any('valence rules give' in x for x in log), log


def test_a_valence_below_the_drawn_sum_is_reported_and_ignored():
    """Two atoms in the corpus state one.  It is not a claim of zero hydrogens, it is nonsense, and
    the arithmetic gives a negative count."""
    _, mol, log = _read([('C', 1), ('C', None), ('C', None)], [(1, 2, 2), (1, 3, 2)])
    assert _counts(mol)[0] == 0
    assert any('cannot be a total valence' in x and 'from the derivation' in x for x in log), log


def test_a_stated_valence_equal_to_the_drawn_sum_is_believed_where_the_rules_are_silent():
    """The one direction the stated valence is taken in, and the case that keeps it in the model:
    "nothing is left over for hydrogen" reads the same under any valence model the writer had."""
    ctab, mol, log = _read([('Fe', 2), ('C', None), ('C', None)], [(1, 2, 1), (1, 3, 1)])
    assert _counts(mol)[0] == 0
    assert not ctab.unknown_hydrogens, log


def test_a_valence_exceeding_the_drawn_sum_where_the_rules_are_silent_is_unknown_not_hydrogen():
    """The 15th rules-silent corpus atom is an iron with valence 3 and one bond drawn; believing the
    difference invents an iron dihydride, so the count is ``H_UNKNOWN`` and the id is in
    ``unknown_hydrogens`` -- a visible refusal rather than something indistinguishable from a bare iron.
    """
    ctab, mol, log = _read([('Fe', 3), ('C', None)], [(1, 2, 1)])
    sid = next(iter(mol.atom_numbers))
    assert sid in ctab.unknown_hydrogens, log
    assert any('coordination or oxidation-state' in x for x in log), log


def test_a_stated_valence_is_not_read_on_an_atom_holding_an_aromatic_bond():
    """The authority order: ``MRV_IMPLICIT_H`` states a count as a count and outranks everything, failing
    that only chython's own chemistry answers, and ``vvv`` is consulted for non-aromatics and nowhere else
    -- scoped to atoms whose bonds all have an integral order, by decision rather than by arithmetic.

    On an aromatic atom the subtraction cannot be repaired here either: ``_drawn_sum`` sees ``{4, 4}`` and
    reaching the real total valence of 3 needs the valence tables and the aromatic classifier, so an
    arithmetic that got it right would be a second copy of the classifier in a format module.  Consulting
    ``vvv`` there also displaces the accurate diagnostic -- a conformant ``vvv=3`` on a pyrrole nitrogen
    reports a valence exceeding the drawn sum "by -5" instead of the ring deciding the atom's class.
    """
    ring = [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 1, 4)]
    bare = [('C', None), ('C', None), ('C', None), ('N', None), ('C', None)]
    stated = [('C', None), ('C', None), ('C', None), ('N', 3), ('C', None)]
    ctab, mol, _ = _read(bare, ring)
    ctab_v, mol_v, log = _read(stated, ring)

    n = [s for s in mol.atom_numbers if mol.element_of(s) == 7][0]
    assert mol.implicit_h_of(n) is None and mol_v.implicit_h_of(n) is None, \
        'a total valence cannot settle the class either way -- see the invariance test below'
    assert ctab.unknown_hydrogens == ctab_v.unknown_hydrogens, 'the field moves no atom'

    # the accurate diagnostic survives the presence of the field, which is the whole fix
    assert any('only the ring decides' in x for x in log), log
    assert not any('exceeds the' in x or 'coordination or oxidation-state' in x for x in log), \
        'the nonsense subtraction is gone, not merely outranked'
    # and the declined field is admitted rather than passed over: the file did say something true
    declined = [x for x in log if 'stated valence 3' in x]
    assert len(declined) == 1, log
    assert str(declined[0]).startswith(UNSUPPORTED), \
        'a legal statement this reader declines to read is our limitation, as `hhh` already is'


def test_the_two_readings_of_an_ambiguous_aromatic_atom_share_one_total_valence():
    """Why no stated valence could settle the class: the two readings trade one unit of ring bond order
    against one hydrogen, so they reach the same total.  Imidazole carries both classes on one ring and
    both nitrogens total 3, so ``vvv`` is invariant across exactly the distinction it would have to
    resolve.  The one case where the arithmetic would discriminate -- the ring-double reading needing
    H = -1 -- the classifier already settles with no ``vvv`` at all, which is the second half below.
    """
    # drawn Kekule, so the classifier is not in the question: N1 donates its lone pair and carries a
    # hydrogen, N3 takes the ring double bond and carries none.  Built rather than parsed from SMILES
    # because nothing in this package may import the facade.
    imidazole = [('N', None), ('C', None), ('N', None), ('C', None), ('C', None)]
    _, m, log = _read(imidazole, [(1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 1, 1)])
    ns = [s for s in m.atom_numbers if m.element_of(s) == 7]
    assert sorted(m.implicit_h_of(s) for s in ns) == [0, 1], f'one class each: {log}'
    totals = {s: sum(m.order_of(s, o) for o in m.neighbors_of(s)) + m.implicit_h_of(s) for s in ns}
    assert set(totals.values()) == {3}, f'and one total valence between them: {totals}'

    # the discriminating case needs no stated valence: N-methylpyrrole's nitrogen has no room for a
    # ring double bond, so the classifier settles it alone and a `vvv` would add nothing.
    ring = [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 1, 4)]
    atoms = [('C', None), ('C', None), ('C', None), ('N', None), ('C', None), ('C', None)]
    ctab, mol, log = _read(atoms, ring + [(4, 6, 1)])
    n = [s for s in mol.atom_numbers if mol.element_of(s) == 7][0]
    assert mol.implicit_h_of(n) == 0 and not ctab.unknown_hydrogens, log


def test_the_reader_declines_on_the_same_predicate_the_corpus_sweep_licenses_itself_with():
    """One predicate, not two spellings of one.  The read gate and the corpus sweep's licence for
    committing a count where chython 2 says unknown have to be the same question, so a change to
    either side's scope moves both and cannot move one silently."""
    from .._hydrogens import _holds_aromatic_bond
    assert holds_an_aromatic_bond is _holds_aromatic_bond


def test_a_stated_valence_on_an_atom_with_no_bonds_drawn_sets_the_count():
    """The one place a stated valence beats a *derived* count.  What disqualifies ``vvv`` elsewhere is that
    it counts orders the file does not always draw; an atom with no bonds has no drawing to be incomplete,
    so ``vvv - drawn == vvv`` and hydrogen is all it can be.  ``AsH3`` as one atom with ``vvv`` 3 is the
    case: the valence collection's free-atom row for every metal and metalloid -- no bonds, no hydrogens,
    a legal state -- answers 0, correctly, and only the reader knows the file said 3.
    """
    ctab, mol, log = _read([('As', 3)], [])
    sid = next(iter(mol.atom_numbers))
    assert _counts(mol) == [3], 'arsine, not a bare arsenic atom'
    assert not ctab.unknown_hydrogens
    assert any('no bonds drawn' in x and 'preference to the 0' in x for x in log), log
    assert not any(str(x).startswith(UNSUPPORTED) for x in log), \
        'the file was fine and the statement was applied; nothing here went unmodelled'


def test_the_stated_valence_of_a_bond_free_atom_is_read_for_every_element():
    """Twenty-one elements, one drawing each: the answer must not depend on whether an element's ladder row
    precedes its free-atom row in the collection, which is an accident of table order.
    """
    for symbol, valence in (('C', 4), ('N', 3), ('O', 2), ('B', 3), ('P', 3), ('S', 2), ('Si', 4),
                            ('Ge', 4), ('Se', 2), ('Te', 2), ('As', 3), ('Sb', 3), ('Bi', 3),
                            ('Al', 3), ('Ga', 3), ('In', 3), ('Tl', 1), ('Sn', 4), ('Pb', 4),
                            ('Po', 2), ('At', 1)):
        ctab, mol, log = _read([(symbol, valence)], [])
        assert _counts(mol) == [valence], f'{symbol} with a stated valence of {valence}: {log}'
        assert not ctab.unknown_hydrogens, symbol


def test_a_bond_free_atom_with_nothing_stated_keeps_its_derived_count():
    """Why the fix is not in the table: making the collection's free-atom row abstain would answer the case
    above at the price of turning every counterion in every SDF into an unknown.  A sodium with nothing
    stated is a sodium ion, so the widening is gated on a statement in the file.
    """
    ctab, mol, log = _read([('Na', None), ('Fe', None)], [])
    assert _counts(mol) == [0, 0]
    assert not ctab.unknown_hydrogens and log == [], log


def test_the_v2000_spelling_of_a_stated_zero_valence_is_read_as_no_hydrogens():
    """``vvv`` 15 is the format's "stated, and it is nothing": the field's own 0 means "not stated", so 15 is
    the only way a V2000 file can say an atom has no bonds and no hydrogens, and it is what this library's
    writer emits for one.  Unread, a bare carbon written as ``vvv`` 15 comes back as methane.
    """
    from .._hydrogens import ZERO_VALENCE

    ctab, mol, log = _read([('C', ZERO_VALENCE)], [])
    assert _counts(mol) == [0], 'a carbon atom, which is what the file said'
    assert any('no bonds drawn' in x and 'preference to the 4' in x for x in log), log


def test_bare_al_and_stated_al_differ_by_one_field_yet_give_opposite_hydrogen_counts():
    """On Al exactly one field separates a metal ion (0 H) from its hydride (3 H): the free-atom row every
    metal and metalloid has — no bonds, no hydrogens, a legal state — reads a bare Al as the ion with no
    log, while a stated valence of 3 on the same atom can only mean three hydrogens.  Both halves are in
    one function so neither can be changed invisibly.

    Two zero cases are here as well.  On Al a stated ``vvv`` 15 agrees with the derivation, so the override
    does not fire and the log cannot distinguish it from the bare case -- the collection's output is an
    integer with no "stated" bit.  Bare carbon with ``vvv`` 15 makes the path observable, its free-atom
    default being 4: losing the stated fact turns the carbon into methane.
    """
    from .._hydrogens import ZERO_VALENCE

    # bare Al, nothing stated: the free-atom row answers 0, no override log
    ctab, mol, log = _read([('Al', None)], [])
    assert _counts(mol) == [0], 'metal: the free-atom row answers 0'
    assert not any('no bonds drawn' in x for x in log), 'nothing was stated, so no override'
    assert not ctab.unknown_hydrogens

    # Al with stated valence 3: the bond-free override selects alane over the free-atom default
    ctab, mol, log = _read([('Al', 3)], [])
    assert _counts(mol) == [3], 'alane: stated valence beats the free-atom 0'
    assert any('no bonds drawn' in x and 'preference to the 0' in x for x in log), log
    assert not ctab.unknown_hydrogens

    # Al with stated-zero (vvv = 15): 0 for a stated reason; indistinguishable from bare in the log
    # because the derivation and the stated value agree — the override branch never fires
    ctab, mol, log = _read([('Al', ZERO_VALENCE)], [])
    assert _counts(mol) == [0], 'stated-zero: the file said no valence, so no hydrogens'
    assert not any('no bonds drawn' in x for x in log), \
        'derivation and stated value both give 0 for Al; neither case emits an override message'
    assert not ctab.unknown_hydrogens

    # C with stated-zero (vvv = 15): 0 even though the free-atom default is 4; the log marks the
    # override, and both the vvv and the MRV_IMPLICIT_H channel survive the round trip
    ctab, mol, log = _read([('C', ZERO_VALENCE)], [])
    assert _counts(mol) == [0], 'stated-zero carbon: the file said no valence, not methane'
    assert any('no bonds drawn' in x and 'preference to the 4' in x for x in log), log
    assert not ctab.unknown_hydrogens
    lines, _ = emit_v2000(mol)
    mol2, _, _ = parse_v2000(lines, []).build()
    assert _counts(mol2) == [0], 'the stated zero survived the round trip'


def test_a_stated_valence_survives_a_round_trip_through_the_bond_free_channel():
    """Written back and read again, arsine is still arsine: both write channels fire for this atom --
    ``vvv`` and the ``MRV_IMPLICIT_H`` group -- and they must say the same thing, so the trip is exact
    whichever one the reader on the other side honours.
    """
    _, mol, _ = _read([('As', 3)], [])
    lines, _ = emit_v2000(mol)
    assert lines[4][48:51].strip() == '3', f'the valence is stated back: {lines[4]!r}'
    mol2, _, _ = parse_v2000(lines, []).build()
    assert _counts(mol2) == [3]


def test_hcount_is_not_in_the_authority_list_at_all():
    """V2000 ``hhh`` and V3000 ``HCOUNT=`` are query fields -- "n or more" -- not statements.  A
    reader that read one as an exact count would disagree with every other tool on the same file."""
    lines = _record([('C', None)], [])
    lines[4] = '    0.0000    0.0000    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0'
    mol, _, log = parse_v2000(lines, []).build()
    assert _counts(mol) == [4]
    assert any('minimum' in x for x in log), log


def test_an_out_of_range_stated_count_is_recomputed_rather_than_stored():
    extra = ['M  STY  1   1 DAT', 'M  SAL   1  1   1', f'M  SDT   1 {MRV_IMPLICIT_H}',
             'M  SED   1 IMPL_H99']
    _, mol, log = _read([('C', None)], [], extra)
    assert _counts(mol) == [4]
    assert any('out of range' in x for x in log), log


def test_a_stated_count_of_fifteen_is_out_of_range_because_fifteen_is_the_sentinel():
    """15 fits the core's H nibble and is still not a count: it is ``H_UNKNOWN``, which is why ``H_MAX`` is
    14.  Split from the ``IMPL_H99`` case above because 99 is refused by any bound while 15 is refused only
    by the right one -- bounded by the nibble's width, a file claiming fifteen hydrogens comes back
    indistinguishable from an atom nobody could compute.
    """
    assert H_MAX == 14, 'a bound that admits 15 admits the sentinel as a count'
    extra = ['M  STY  1   1 DAT', 'M  SAL   1  1   1', f'M  SDT   1 {MRV_IMPLICIT_H}',
             f'M  SED   1 IMPL_H{H_MAX + 1}']
    ctab, mol, log = _read([('C', None)], [], extra)
    assert _counts(mol) == [4], 'recomputed from the rules, not believed'
    assert any('out of range' in x for x in log), log
    assert not ctab.unknown_hydrogens, \
        'the rules answered for this carbon, so a rejected statement does not make it unknown'


def test_no_path_through_calc_implicit_can_store_the_sentinel_value(corpus):
    """The invariant behind :data:`H_MAX`, over every atom of the corpus: three paths write a count -- a
    stated one, the valence rules, a stated valence where the rules are silent -- and one of them admitting
    15 makes the sentinel reachable.
    """
    seen = 0
    for records in corpus.values():
        for record in records:
            mol = parse_record(record)
            for sid in mol.atom_numbers:
                h = mol.implicit_h_of(sid)
                assert h is None or 0 <= h <= H_MAX, f'atom {sid} stored {h}'
                seen += 1
    assert seen > 5000, seen


# ----------------------------------------------------------------------------------- the unknowns

def test_an_atom_no_rule_covers_is_unknown_and_not_a_refusal():
    """A five-bonded carbon.  Real editors write these, and a reader does not get to reject one:
    permissive in, permissive to store, honest when asked."""
    ctab, mol, log = _read([('C', None)] + [('F', None)] * 5,
                           [(1, i, 1) for i in range(2, 7)])
    sid = next(iter(mol.atom_numbers))
    assert sid in ctab.unknown_hydrogens
    # the one message in the module that earns `unsupported: `: the valence collection describes nothing
    # for a five-bonded carbon, so the limitation is chython's and the file may be well formed.  Its
    # counterpart, an aromatic atom only the ring can settle, carries no prefix.
    assert any(str(x).startswith(f'{UNSUPPORTED}atom {sid}:') and 'no valence rule' in x
               and 'not known' in x for x in log), log
    assert any('unknown implicit hydrogen count' in x for x in log), \
        'the summary line is what makes the question answerable by one grep of the log'


def test_every_count_the_reader_produces_is_assigned_and_never_inherited():
    """The reader's counts are a function of the *record*, never of ``add_atom``'s default: ``calc_implicit``
    walks ``atom_numbers`` and every branch writes a count before it continues.

    Pinned by making the default observable rather than by naming its value, so the test needs no edit when
    the default changes and cannot pass by agreeing with it.  Two fixtures, because neither covers both
    defaults: methane's carbon is 4, so an inherited 0 or 15 both show; carbon tetrafluoride's carbon is a
    genuine 0, where only a sentinel default shows.  The ``!= default`` assertion runs first so a leak
    reports which hazard it hit.
    """
    from ....core import MoleculeContainer

    bare = MoleculeContainer()
    default = bare.implicit_h_of(bare.add_atom(6))

    ctab, mol, _ = _read([('C', None)], [])
    assert _counts(mol) != [default], f'inherited the add_atom default ({default}) instead of asking'
    assert _counts(mol) == [4], 'a lone carbon in a file is methane'

    ctab, mol, _ = _read([('C', None)] + [('F', None)] * 4, [(1, i, 1) for i in range(2, 6)])
    assert _counts(mol) == [0, 0, 0, 0, 0], 'carbon tetrafluoride: a real zero, on every atom'
    assert not ctab.unknown_hydrogens, 'and determinable, so nothing here is the sentinel'


def test_an_unknown_count_does_not_leave_the_arena_as_fifteen_hydrogens():
    """An unknown count must not reach a derived representation *as a number*: a consumer reading the nibble
    raw turns "nobody knows" into "there are fifteen" silently, as ``molecule_to_inchi`` once did -- alanine
    with one unresolved carbon came back ``C3H19NO2``, ``1H15`` in the h layer, a wrong formula under a
    wrong InChIKey.  ``_inchi.pxi`` writes ``-1``, InChI's own "auto", and nothing but this asserts the
    sentinel is not handed to libinchi as a number.
    """
    from ....core import inchi_library_loaded, molecule_to_inchi

    if not inchi_library_loaded():
        skip('no libinchi for this platform')

    # A five-bonded carbon carrying five fluorines: the rules have nothing to say, so its count is
    # unknown, and the rest of the molecule is perfectly ordinary.
    ctab, mol, _ = _read([('C', None)] + [('F', None)] * 5, [(1, i, 1) for i in range(2, 7)])
    assert len(ctab.unknown_hydrogens) == 1, 'the fixture has to actually produce one'
    out = molecule_to_inchi(mol)
    layers = [x for x in out.split('/') if x.startswith('h')]
    assert 'H15' not in out and not any('H15' in x for x in layers), \
        f'the sentinel went out as a hydrogen count: {out}'


def test_the_unknown_set_is_one_attribute_and_not_a_walk_over_atoms():
    """The repair pipeline's first question is "does this record contain counts nobody stated", asked
    once per record over a large SDF.  It has to be answerable without touching every atom."""
    ctab, mol, _ = _read([('C', None), ('O', None)], [(1, 2, 1)])
    assert ctab.unknown_hydrogens == ()
    ctab, mol, _ = _read([('C', None)] + [('F', None)] * 5, [(1, i, 1) for i in range(2, 7)])
    assert len(ctab.unknown_hydrogens) == 1


def test_where_the_count_is_unknown_the_writer_stays_silent():
    """Writing ``IMPL_H0`` for an unknown count turns "nobody said" into "there are none", asserting
    in the file the very fact the reader was careful not to invent.  A silent file stays silent."""
    ctab, mol, _ = _read([('C', None)] + [('F', None)] * 5, [(1, i, 1) for i in range(2, 7)])
    lines, log = emit_v2000(mol)
    assert not any(MRV_IMPLICIT_H in x for x in lines), lines
    assert lines[4][48:51].strip() in ('', '0'), f'no valence stated either: {lines[4]!r}'
    assert any('no known implicit hydrogen count' in x for x in log), log


def test_the_writer_stays_silent_for_a_molecule_it_never_read_from_a_file():
    """This molecule did not come through the reader, so no caller could hand the writer a set of unknown
    atoms.  It is silent anyway because the fact is in the storage: ``implicit_h_of`` answers ``None`` and
    :func:`valence_for_write` follows.  A receipt travelling beside the molecule instead of inside it is
    what made a hand-built or transformed molecule go out claiming ``IMPL_H0``.
    """
    from ....core import H_UNKNOWN, MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        sid = mol.add_atom('C')
        mol.set_hydrogens(sid, H_UNKNOWN)
    assert mol.implicit_h_of(sid) is None
    assert valence_for_write(mol, sid) is None, 'nothing to state, so nothing is stated'

    lines, log = emit_v2000(mol)
    assert not any(MRV_IMPLICIT_H in x for x in lines), lines
    assert lines[4][48:51].strip() in ('', '0'), f'no valence stated either: {lines[4]!r}'
    assert any('no known implicit hydrogen count' in x for x in log), \
        'the loss is still reported; a writer that drops a fact silently is how it disappears'


# ------------------------------------------------------------------------ the round-trip channel

def test_a_count_the_rules_disagree_with_survives_a_round_trip():
    """The reason ``MRV_IMPLICIT_H`` is written at all.  ``vvv`` cannot carry a count back into this
    library, because this library ranks it below the rules and would override it on the way in."""
    ctab, mol, _ = _read([('C', None), ('O', None)], [(1, 2, 1)])
    sid = next(iter(mol.atom_numbers))
    with mol.edit():
        mol.set_hydrogens(sid, 1)  # not what the rules give
    lines, log = emit_v2000(mol)
    assert any(MRV_IMPLICIT_H in x for x in lines), lines
    mol2, _, _ = parse_v2000(lines, []).build()
    assert mol2.implicit_h_of(next(iter(mol2.atom_numbers))) == 1


def test_the_writer_is_idempotent_and_does_not_accumulate_records():
    """Read, write, read, write must not add a data S-group per pass.  Records are re-derived from
    the stored counts rather than passed through, which is what makes that true."""
    ctab, mol, _ = _read([('C', None), ('O', None)], [(1, 2, 1)])
    sid = next(iter(mol.atom_numbers))
    with mol.edit():
        mol.set_hydrogens(sid, 1)
    lines, _ = emit_v2000(mol)
    for _ in range(3):
        ctab2 = parse_v2000(lines, [])
        mol2, store2, _ = ctab2.build()
        lines, _ = emit_v2000(mol2, store2)
    assert sum(1 for x in lines if MRV_IMPLICIT_H in x) == 1, lines


def test_valence_for_write_is_silent_when_the_rules_already_give_the_count():
    """A record stating a valence on every atom is noise, and noise other tools then honour."""
    _, mol, _ = _read([('C', None), ('O', None)], [(1, 2, 1)])
    assert all(valence_for_write(mol, s) is None for s in mol.atom_numbers)


def test_a_genuine_zero_valence_is_written_as_the_format_s_own_spelling():
    """The field's 0 means "not stated", so a valence that really is nothing needs 15."""
    from .._hydrogens import ZERO_VALENCE
    from ....core import MoleculeContainer

    mol = MoleculeContainer()
    with mol.edit():
        sid = mol.add_atom('C')
        mol.set_hydrogens(sid, 0)
    assert valence_for_write(mol, sid) == ZERO_VALENCE


# ------------------------------------------------------------------------------- the unit helpers

def test_implicit_from_valence_returns_none_for_a_statement_that_cannot_be_one():
    _, mol, _ = _read([('C', None), ('C', None)], [(1, 2, 2)])
    sid = next(iter(mol.atom_numbers))
    assert implicit_from_valence(mol, sid, 4) == 2
    assert implicit_from_valence(mol, sid, 1) is None, 'below the drawn sum'
    assert implicit_from_valence(mol, sid, H_MAX + 3) is None, 'more hydrogens than an atom can hold'


def test_implicit_for_atom_distinguishes_no_rule_from_a_rule_saying_zero():
    """No rule against a rule saying zero -- the distinction the third state rests on.  Both halves are
    asserted on the reason code as well as the count, since ``(None, HYD_NO_VALENCE_RULE)`` and
    ``(None, HYD_AMBIGUOUS_AROMATIC)`` are the same non-answer for opposite reasons and the log prefix
    hangs off which one it is.
    """
    _, mol, _ = _read([('C', None)] + [('F', None)] * 4, [(1, i, 1) for i in range(2, 6)])
    sids = list(mol.atom_numbers)
    assert implicit_for_atom(mol, sids[0]) == (0, HYD_DERIVED), \
        'four bonds: a rule, and it says zero'
    _, mol, _ = _read([('C', None)] + [('F', None)] * 5, [(1, i, 1) for i in range(2, 7)])
    sids = list(mol.atom_numbers)
    assert implicit_for_atom(mol, sids[0]) == (None, HYD_NO_VALENCE_RULE), \
        'five bonds: no rule at all'


def test_implicit_for_atom_passes_count_stated_false_and_that_is_the_whole_ruling():
    """The one thing this module adds to the core's derivation: a CTfile bond block entry of type 4 has no
    channel for pyrrole-versus-pyridine, so that pnictogen comes back ``None`` here while the same atom
    asked with ``count_stated=True``, which is what SMILES passes, answers 0.  Both calls, side by side.
    """
    from ....core._core import derive_implicit_hydrogen

    # pyridine drawn all-aromatic: the nitrogen's class is the ring's to decide
    _, mol, _ = _read([('N', None)] + [('C', None)] * 5,
                      [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4)])
    nitrogen = next(iter(mol.atom_numbers))
    assert implicit_for_atom(mol, nitrogen) == (None, HYD_AMBIGUOUS_AROMATIC)
    assert derive_implicit_hydrogen(mol, nitrogen, count_stated=True) == (0, HYD_DERIVED), \
        'the argument is what makes the difference, so it has to be observable through it'


def test_hydrogen_carries_no_hydrogen():
    _, mol, _ = _read([('H', None), ('C', None)], [(1, 2, 1)])
    assert _counts(mol)[0] == 0


def test_hydrogen_carries_no_hydrogen_however_it_is_drawn():
    """Four drawings the valence collection has no row for, and hydrogen answers 0 for all of them.  The
    collection describes hydrogen at a drawn order sum of one, the only sum a well drawn hydrogen has;
    files carry the others anyway -- a lone atom, a single dative contact, a diborane bridge -- and the
    lookup returns ``None`` for each.  Zero is not a fallback but a property of the element, so the core
    states it ahead of the lookup and this reader does not state it again.

    Six of the corpus's 438 explicit hydrogens are drawn at a sum other than one, so losing that element
    fact trades six correct zeros for six sentinels.  Each drawing is asserted separately, since they fail
    for different arithmetic.
    """
    for label, atoms, bonds in (
            ('a lone atom', [('H', None)], []),
            ('one dative contact', [('H', None), ('B', None)], [(1, 2, 8)]),
            ('a diborane-style bridge', [('H', None), ('B', None), ('B', None)],
             [(1, 2, 1), (1, 3, 1)]),
            ('two dative contacts', [('H', None), ('B', None), ('B', None)],
             [(1, 2, 8), (1, 3, 8)])):
        ctab, mol, _ = _read(atoms, bonds)
        sid = next(iter(mol.atom_numbers))
        assert implicit_for_atom(mol, sid) == (0, HYD_DERIVED), label
        assert mol.implicit_h_of(sid) == 0, label
        assert not ctab.unknown_hydrogens, f'{label}: a hydrogen turned into a sentinel'


def test_every_explicit_hydrogen_in_the_corpus_is_answered_and_six_of_them_need_no_table(corpus):
    """The population behind the test above: 438 explicit hydrogens, of which 6 are drawn at a bond order
    sum the collection has no row for -- 4 at a sum of 0, 2 at a sum of 2.  Both halves are asserted, every
    hydrogen getting a 0 and the six existing, so the answer is not coming from a row.
    """
    from ....core._core import derive_implicit_hydrogen

    total = off_the_table = 0
    for records in corpus.values():
        for record in records:
            mol = parse_record(record)
            for sid in mol.atom_numbers:
                if mol.element_of(sid) != 1:
                    continue
                total += 1
                assert mol.implicit_h_of(sid) == 0, f'atom {sid} of a corpus record'
                drawn = sum(mol.order_of(sid, o) for o in mol.neighbors_of(sid)
                            if mol.order_of(sid, o) != 8)
                if drawn != 1:
                    off_the_table += 1
                    # asked of the core directly, so the claim is about the shared derivation and not
                    # about this reader's ranking of it
                    assert derive_implicit_hydrogen(mol, sid) == (0, HYD_DERIVED), sid
    assert (total, off_the_table) == (438, 6), (total, off_the_table)


# ------------------------------------------- what the core reports back, and what reaches the log

# One test per reason code, asserting on the *line* and not only on the count, since the count cannot
# distinguish "no rule" from "the ring decides" and the prefix is the whole difference.
# `HYD_NO_VALENCE_RULE`'s test is `test_an_atom_no_rule_covers_is_unknown_and_not_a_refusal` above, where
# the rest of that outcome's behaviour lives.

def test_a_derived_count_reports_nothing():
    """``HYD_DERIVED`` with no flag is the ordinary case and the log stays empty.  Asserted because
    every other reason code here is recognised by the line it adds, and a reader that logged on the
    ordinary case would drown all four of them."""
    _, mol, log = _read([('C', None), ('O', None)], [(1, 2, 1)])
    assert [implicit_for_atom(mol, s)[1] for s in mol.atom_numbers] == [HYD_DERIVED, HYD_DERIVED]
    assert log == [], log


def test_an_aromatic_atom_only_the_ring_can_settle_is_reported_without_the_prefix():
    """``HYD_AMBIGUOUS_AROMATIC``, and the prefix easiest to get backwards: chython models an aromatic
    pyrrole nitrogen perfectly well and what is missing is a statement the CTfile has no channel for, so no
    ``unsupported: ``.  A line rather than silence, because only it says what repair -- ``kekule()``, or an
    ``MRV_IMPLICIT_H`` statement -- answers the question.
    """
    ctab, mol, log = _read([('N', None)] + [('C', None)] * 5,
                           [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4)])
    sid = next(iter(mol.atom_numbers))
    assert implicit_for_atom(mol, sid) == (None, HYD_AMBIGUOUS_AROMATIC)
    assert sid in ctab.unknown_hydrogens and mol.implicit_h_of(sid) is None
    lines = [x for x in log if str(x).startswith(f'atom {sid}:') and 'not known' in x]
    assert lines, log
    assert not str(lines[0]).startswith(UNSUPPORTED), \
        'an ambiguity the format cannot express is the file falling short, not chython'
    assert 'aromatic bond(s)' in lines[0] and 'no Kekule form' in lines[0], lines[0]


def test_a_count_from_the_other_aromatic_reading_is_stored_and_reported_without_the_prefix():
    """``HYD_DERIVED_OTHER_READING``: a neutral phosphorus drawn with two aromatic bonds and an exocyclic
    double bond, where the class the classifier picks has no valence row and the other class does.  There
    is a count, so no ``unsupported: `` -- and the line is not optional either, since ``kekule()`` will
    pick the classifier's class and the stored count will then disagree with the Kekule structure.
    """
    ctab, mol, log = _read([('P', None), ('C', None), ('C', None), ('O', None)],
                           [(1, 2, 4), (1, 3, 4), (1, 4, 2)])
    sid = next(iter(mol.atom_numbers))
    assert implicit_for_atom(mol, sid) == (0, HYD_DERIVED_OTHER_READING)
    assert mol.implicit_h_of(sid) == 0 and sid not in ctab.unknown_hydrogens
    lines = [x for x in log if str(x).startswith(f'atom {sid}:')]
    assert lines and 'other reading' in lines[0] and 'kekule()' in lines[0], log
    assert not any(str(x).startswith(UNSUPPORTED) for x in log), \
        'the count was derived and stored; nothing here went unmodelled'


def test_an_atom_with_no_aromatic_form_is_reported_beside_its_count():
    """``HYD_NO_AROMATIC_FORM`` is a flag OR-ed onto an outcome, not a fifth outcome: a neutral beryllium
    drawn with two aromatic bonds has no aromatic form, is read as saturated, and still gets a count, so the
    observation is reported alongside ``HYD_DERIVED`` and the mask keeps the outcome recognisable.  No
    prefix, the atom being stored exactly as drawn.  The message names the number of aromatic bonds because
    the flag means "this element in this state" -- a carbon with one aromatic bond reaches the same line.
    """
    ctab, mol, log = _read([('Be', None), ('C', None), ('C', None)], [(1, 2, 4), (1, 3, 4)])
    sid = next(iter(mol.atom_numbers))
    h, reason = implicit_for_atom(mol, sid)
    assert (h, reason & HYD_REASON_MASK) == (0, HYD_DERIVED)
    assert reason & HYD_NO_AROMATIC_FORM
    assert mol.implicit_h_of(sid) == 0 and not ctab.unknown_hydrogens
    lines = [x for x in log if str(x).startswith(f'atom {sid}:')]
    assert lines and 'no aromatic form with 2 aromatic bond(s)' in lines[0], log
    assert not any(str(x).startswith(UNSUPPORTED) for x in log), log
    # the same line for the carbons, whose aromatic form exists but not with one bond drawn
    assert sum(1 for x in log if 'no aromatic form with 1 aromatic bond(s)' in x) == 2, log


def test_a_stated_count_does_not_silence_the_observation_about_the_bonds_drawn():
    """A stated count and an unresolvable aromatic system are two facts and the reader owes both:
    ``MRV_IMPLICIT_H`` is rank 1 and describes the *count*, while ``HYD_NO_AROMATIC_FORM`` describes the
    *drawing*, whose bonds are stored as order 4 regardless.  So the top-authority channel must not take
    the observation down with it.  No corpus record does this -- the two flagged corpus atoms state nothing
    -- so only a fixture can hold it.
    """
    extra = ['M  STY  1   1 DAT', 'M  SAL   1  1   1', f'M  SDT   1 {MRV_IMPLICIT_H}',
             'M  SED   1 IMPL_H0']
    atoms, bonds = [('Be', None), ('C', None), ('C', None)], [(1, 2, 4), (1, 3, 4)]
    ctab, mol, silent = _read(atoms, bonds)
    sid = next(iter(mol.atom_numbers))
    assert any(str(x).startswith(f'atom {sid}:') and 'no aromatic form with 2' in x for x in silent), silent

    ctab, mol, log = _read(atoms, bonds, extra)
    assert _counts(mol) == [0, 3, 3], 'the stated count is honoured, which is the whole of rank 1'
    assert mol.aromatic_bond_count == 2, 'and the bonds it says nothing about are stored as drawn'
    assert any(str(x).startswith(f'atom {sid}:') and 'no aromatic form with 2' in x for x in log), log
    assert not any(str(x).startswith(UNSUPPORTED) for x in log), \
        'an observation about the input is not a construct we failed to model'


def test_exactly_one_silent_outcome_can_reach_the_stated_valence_line_and_it_is_always_ours():
    """The prefix on the "valence exceeds the drawn sum" line is a constant, not a run-time decision:
    scoping the valence channel to non-aromatic atoms makes ``HYD_AMBIGUOUS_AROMATIC`` unreachable there by
    construction -- it requires an aromatic bond, and such an atom has its stated valence declined before
    the comparison -- so the line is always ``HYD_NO_VALENCE_RULE``, chython's own gap, and prefixed.  Both
    fixtures stay: a regression shows up as a resurrected line on the pyridine nitrogen.
    """
    # no valence rule (5 drawn single bonds on a carbon), no aromatic bond, valence 6 stated
    ctab, mol, log = _read([('C', 6)] + [('F', None)] * 5, [(1, i, 1) for i in range(2, 7)])
    sid = next(iter(mol.atom_numbers))
    _, reason = implicit_for_atom(mol, sid)
    assert (reason & HYD_REASON_MASK, bool(reason & HYD_NO_AROMATIC_FORM)) \
        == (HYD_NO_VALENCE_RULE, False), reason
    assert sid in ctab.unknown_hydrogens
    lines = [x for x in log if 'no derivable hydrogen count' in x]
    assert lines and str(lines[0]).startswith(f'{UNSUPPORTED}atom {sid}:'), log
    assert 'exceeds the 5 drawn by 1' in lines[0], log

    # the ambiguous-aromatic counterpart, with the same stated valence: the sentence is not reached at
    # all now, and what the atom is told instead is the accurate thing about its ring
    ctab, mol, log = _read([('N', 9)] + [('C', None)] * 5,
                           [(1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 1, 4)])
    sid = next(iter(mol.atom_numbers))
    assert implicit_for_atom(mol, sid) == (None, HYD_AMBIGUOUS_AROMATIC)
    assert sid in ctab.unknown_hydrogens
    assert not [x for x in log if 'no derivable hydrogen count' in x], log
    assert any(str(x).startswith(f'atom {sid}:') and 'only the ring decides' in x for x in log), log
    assert any(str(x).startswith(UNSUPPORTED) and f'stated valence 9 not read' in x for x in log), log


def test_a_reason_code_this_reader_has_no_ruling_for_is_named_and_not_called_unsupported(monkeypatch):
    """The prefix decision defaults towards claiming nothing: the reader recognises four codes and prefixes
    exactly one, so a fifth falls through both arms and is *named* rather than classified -- the value of
    ``unsupported: `` being that a caller screening on it can trust what it means.  Unreachable with the
    core's four codes, ``h is None`` implying one of two, so the derivation is substituted to reach it.
    """
    from .. import _hydrogens

    monkeypatch.setattr(_hydrogens, 'implicit_for_atom', lambda mol, sid: (None, 4))
    _, mol, _ = _read([('C', None)], [])
    result = calc_implicit(mol)
    sid = next(iter(mol.atom_numbers))
    assert mol.implicit_h_of(sid) is None and result.unknown == (sid,), 'still stored, still marked'
    lines = [x for x in result.log if str(x).startswith(f'atom {sid}:')]
    assert lines and 'reason 4' in lines[0] and 'no ruling' in lines[0], result.log
    assert not any(str(x).startswith(UNSUPPORTED) for x in result.log), \
        'an unrecognised code says nothing about whose gap it is, so it may not claim ours'


def test_calc_implicit_returns_a_result_object_that_cannot_be_unpacked():
    """It has already grown from one field to two and will grow again.  A positional read has to fail
    now rather than silently misread later."""
    _, mol, _ = _read([('C', None)], [])
    result = calc_implicit(mol)
    assert hasattr(result, 'log') and hasattr(result, 'unknown')
    try:
        a, b = result
    except TypeError:
        pass
    else:
        raise AssertionError('HydrogenResult became unpackable; a caller can now read it positionally')


def test_every_corpus_record_reads_and_the_unknowns_are_exactly_v2s(corpus, v2_molecules):
    """The whole 512-record corpus reads, and the marking is checked against chython 2 per atom rather than
    by a count:

    * every atom V2 answers ``None`` for is in ``unknown_hydrogens`` -- nothing gets a fabricated count;
    * every atom V2 commits a number for is **not** in ``unknown_hydrogens``, and the number agrees -- which
      is the direction a lazy implementation fails, marking everything unknown satisfying the first.

    The first is deliberately asymmetric with the second.  V2 answers ``None`` for 253 atoms here; this
    reader marks 124 and answers the other 129, because the shared derivation asks the kekuliser's own
    aromatic classifier where V2 asked a narrower arithmetic, so a pyridine-class heteroatom, a charged
    aromatic and a non-organic-subset element get a count every Kekule form of the ring agrees on.
    Recovering a count V2 leaves unknown is allowed; committing a *different* count, or losing one V2
    commits, is not.

    So ``recovered`` is checked for its two legitimate sources -- a stated valence equal to the drawn sum,
    which V2's ``parse_mol_v2000`` never reads, and an aromatic bond.  An atom in neither would mean the two
    valence collections had drifted.
    """
    from .._sdf import sniff_version
    from .._v3000 import V3000_STAMP, parse_v3000

    total = aromatic = no_rule = 0
    v2_unknown_and_marked = v2_number_and_agreed = v2_number_but_marked = disagreed = 0
    recovered = []
    for name, records in corpus.items():
        for n, record in enumerate(records):
            total += 1
            parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
            ctab = parse(record, [])
            mol, _, log = ctab.build()  # no flags: reading must need no opt-in
            unknown = set(ctab.unknown_hydrogens)
            if unknown:
                assert any('unknown implicit hydrogen count' in x for x in log), \
                    'a marked atom without its receipt in the log is what the ruling forbids'
            # every marked atom must name its own reason, not just appear in the summary.  The
            # `unsupported: ` prefix is stripped before the atom id is read off, since whether a line
            # carries one is what distinguishes chython's gap from the file not having said.
            reason = {}
            for line in log:
                line_s = str(line)
                plain = line_s[len(UNSUPPORTED):] if line_s.startswith(UNSUPPORTED) else line_s
                if plain.startswith('atom ') and 'not known' in plain:
                    ring_said = 'aromatic bond(s)' in plain
                    assert ring_said != line_s.startswith(UNSUPPORTED), \
                        f'the prefix and the reason disagree on this line: {line!r}'
                    reason[int(plain.split()[1].rstrip(':'))] = 'aromatic' if ring_said else 'no_rule'
            assert unknown <= reason.keys(), \
                f'marked unknown with no per-atom log line: {unknown - reason.keys()}'
            aromatic += sum(1 for s in unknown if reason[s] == 'aromatic')
            no_rule += sum(1 for s in unknown if reason[s] == 'no_rule')
            if name not in v2_molecules:
                continue
            m2 = v2_molecules[name][n]
            # V2 numbers its atoms 1..n in atom-block order and so do our stable ids on a fresh build
            for i, (sid, num) in enumerate(zip(mol.atom_numbers, m2)):
                h2 = m2.atom(num).implicit_hydrogens
                if h2 is None:
                    if sid in unknown:
                        v2_unknown_and_marked += 1
                    else:
                        recovered.append((name, n, i, ctab.atoms[i].valence,
                                          holds_an_aromatic_bond(mol, sid)))
                elif sid in unknown:
                    v2_number_but_marked += 1
                elif mol.implicit_h_of(sid) == h2:
                    v2_number_and_agreed += 1
                else:
                    disagreed += 1
    assert total == 512, f'the corpus changed size ({total}); remeasure before trusting the rest'
    assert v2_number_but_marked == 0, \
        f'{v2_number_but_marked} atom(s) marked unknown that V2 answers outright'
    assert disagreed == 0, f'{disagreed} atom(s) where both commit and the numbers differ'
    assert (aromatic, no_rule) == (37, 87), (aromatic, no_rule)
    assert v2_unknown_and_marked == aromatic + no_rule, \
        'every atom marked unknown here is one V2 answers None for; the converse no longer holds'

    # Two ways to earn a recovery, and the split is asserted rather than the total so a recovery that is
    # neither shows up as drift.  Fourteen come from a stated total valence equal to the drawn bond sum,
    # which V2's `parse_mol_v2000` never reads (`line[48:51]`); 129 come from an aromatic bond, where the
    # shared derivation asks the kekuliser's classifier and V2 asked its own narrower arithmetic.
    assert len(recovered) == 143, len(recovered)
    assert sum(1 for *_, valence, _ in recovered if valence is not None) == 14
    assert all(valence is not None or aromatic for *_, valence, aromatic in recovered), \
        f'a count committed where V2 says unknown, with neither a stated valence nor an aromatic bond'


def test_implicit_h_records_skips_the_atoms_that_hold_the_sentinel():
    """Same two cases as before, with the fact moved from an argument into the atom.  1 is a count and
    gets a record; ``H_UNKNOWN`` is not a count and gets nothing."""
    from ....core import H_UNKNOWN

    _, mol, _ = _read([('C', None), ('O', None)], [(1, 2, 1)])
    sid = next(iter(mol.atom_numbers))
    with mol.edit():
        mol.set_hydrogens(sid, 1)
    records, _ = implicit_h_records(mol)
    assert [r.name for r in records] == [MRV_IMPLICIT_H]
    with mol.edit():
        mol.set_hydrogens(sid, H_UNKNOWN)
    records, _ = implicit_h_records(mol)
    assert records == []
