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
"""The SMILES reader.

The invariants these tests exist to protect, in order of how expensive they are to lose:

  * aromatic input is stored aromatic and the reader NEVER kekulises
  * a Kekule input is stored exactly as written
  * a syntax error raises with an offset; chemistry never raises
  * atom-case aromatic promotion fires on an all-lowercase ring and NOT on biphenyl's inter-ring
    bond, and every promotion is logged
  * `@`, `/` and `\\` mean what the string said: the sign is measured against RDKit's InChI, and a
    configuration the reader cannot place is NAMED in the log rather than dropped
"""
from pytest import mark, raises

from chython.core import MoleculeContainer
from chython.core._core import (IncorrectSmiles, STEREO_ABS, STEREO_AND, STEREO_OR, kekule,
                                read_smiles, smi_element_table, write_smiles)

try:
    from chython.core import inchi_library_loaded, molecule_to_inchi
    _INCHI = inchi_library_loaded()
except ImportError:                                  # a build that predates `_inchi.pxi`
    _INCHI = False
    molecule_to_inchi = None

# The absolute sign of a configuration cannot be checked without something outside this tree that
# already knows it, so those tests -- and only those -- skip when the oracle is absent.
needs_inchi = mark.skipif(not _INCHI, reason='libinchi not loaded, so the external oracle is absent')


def orders(mol):
    """{(n, m): order} with n < m, so a test can state bonds without caring about traversal."""
    out = {}
    for b in mol.bonds():
        out[(b.n, b.m) if b.n < b.m else (b.m, b.n)] = b.order
    return out


def hydrogens(mol):
    out = []
    for a in mol.atoms():
        out.append(a.implicit_h)
    return out


def elements(mol):
    out = []
    for a in mol.atoms():
        out.append(a.element)
    return out


# ---------------------------------------------------------------- the element table


def test_element_table_is_the_same_authority_as_the_symbol_table():
    """SMI_ELEMENT against `add_atom`, which reads SYMBOL_TO_NUMBER, which is built from SYMBOLS.

    Not a spot check: every entry of the generated table is resolved through the authority, and the
    values are required to cover 1..118 exactly.  A table that drifted from `_elements.pxi` -- an
    element renamed, a row shifted -- fails here rather than in a parse three months later.
    """
    table = smi_element_table()
    assert len(table) == 118
    assert sorted(table.values()) == list(range(1, 119))
    mol = MoleculeContainer()
    with mol.edit():
        for symbol, number in table.items():
            mol.add_atom(symbol)
    assert elements(mol) == list(table.values())


def test_two_letter_symbols_the_chython_two_regex_cannot_read():
    # chython 2's atom_re excludes `c` from a symbol's second position, so `[Sc]` is unreadable
    # there.  It is scandium here, and greedy matching inside a bracket is why.
    assert elements(read_smiles('[Sc]')) == [21]
    assert elements(read_smiles('[Se]')) == [34]
    assert elements(read_smiles('[Hg]')) == [80]
    assert elements(read_smiles('[Og]')) == [118]


def test_outside_a_bracket_the_match_is_not_greedy():
    # `SC` is sulfur bonded to carbon and must never read as scandium: outside brackets only the
    # organic subset exists, and only Cl and Br are two characters long
    assert elements(read_smiles('SC')) == [16, 6]
    assert elements(read_smiles('ClCBr')) == [17, 6, 35]
    assert elements(read_smiles('BC')) == [5, 6]


# ---------------------------------------------------------------- bracket atoms


def test_bracket_fields():
    mol = read_smiles('[13CH3-]')
    atom = next(mol.atoms())
    assert (atom.element, atom.isotope, atom.implicit_h, atom.charge) == (6, 13, 3, -1)
    mol = read_smiles('[Fe+2]')
    assert next(mol.atoms()).charge == 2
    # the repeated-sign spelling means the same thing
    assert next(read_smiles('[Fe++]').atoms()).charge == 2
    assert next(read_smiles('[O-]').atoms()).charge == -1
    assert next(read_smiles('[NH4+]').atoms()).implicit_h == 4
    assert next(read_smiles('[C:12]').atoms()).map_number == 12


def test_bracket_field_order_is_not_enforced():
    # OpenSMILES fixes the order; writers in the wild do not agree with it, and a reader that
    # refuses `[C-H3]` refuses input it could have understood
    for text in ('[CH3-]', '[C-H3]'):
        atom = next(read_smiles(text).atoms())
        assert (atom.implicit_h, atom.charge) == (3, -1)


def test_a_bracket_always_states_its_hydrogen_count():
    # `[CH0]` and `[C]` are the same statement -- zero -- and both differ from bare `C`, which
    # states nothing and gets the notation model's four
    assert hydrogens(read_smiles('[C]')) == [0]
    assert hydrogens(read_smiles('[CH0]')) == [0]
    assert hydrogens(read_smiles('C')) == [4]


def test_a_bracket_label_is_the_marker_carrying_that_text():
    # `[Pol]`, `[Resin]`, `[REG42]`: a display label, a polymer end, a registry identifier. Element 0
    # plus the text as the atom's alias -- and `[Pol]` is the case that matters, because `Po` IS an
    # element and the `l` after it is what says this bracket names something else.
    for text, label in (('[Pol]CC', b'Pol'), ('[Resin]C', b'Resin'), ('[REG42]', b'REG42'),
                        ('[OMe]C', b'OMe'), ('[h11b6]C', b'h11b6'), ('[Zz]', b'Zz')):
        log = []
        mol = read_smiles(text, log)
        atom = next(a for a in mol.atoms() if a.is_r)
        assert atom.element == 0 and atom.r_index == 0, text
        assert mol.aliases == {atom.n: label}, text
        assert len(log) == 1 and log[0].rule == 'smiles:label-as-marker', (text, log)
        assert log[0].severity == 'info', text
    # the bracket's other fields still read, since the label run stops at every one of them
    mol = read_smiles('[Rgp+:7]C')
    atom = next(a for a in mol.atoms() if a.is_r)
    assert (atom.charge, atom.map_number) == (1, 7)
    assert mol.aliases == {atom.n: b'Rgp'}
    # and an element followed by a field is untouched
    for text in ('[NH4+]', '[13CH3-]', '[Po]', '[Se]', '[R12]C', '[nH]1cccc1'):
        log = []
        mol = read_smiles(text, log)
        assert not mol.aliases, text
        assert not [r for r in log if r.rule == 'smiles:label-as-marker'], text


def test_a_charge_outside_the_field_is_clamped_and_reported_lost():
    # `[Pt+10]`, `[ZrH8+12]`: a writer that turned every dative contact into a formal charge pair. The
    # connectivity is worth having, so the charge is clamped to the field and the loss is named.
    log = []
    mol = read_smiles('[Pt+10]', log)
    assert next(mol.atoms()).charge == 8
    assert len(log) == 1 and 'charge 10' in log[0] and 'stored as 8' in log[0]
    log = []
    assert next(read_smiles('[O-5]', log).atoms()).charge == -4
    assert len(log) == 1 and 'charge -5' in log[0]
    # the repeated-sign spelling clamps the same way, and inside the field nothing is logged
    log = []
    assert next(read_smiles('[O--]', log).atoms()).charge == -2
    assert not log


def test_a_doubled_directional_token_is_one_token():
    # `C\\3`, `)\\N2`: an unescape the producer lost, and only the backslash doubles. The pair reads
    # as the direction it stood for; a third in a row is still two bond tokens in a row.
    log = []
    mol = read_smiles(r'C/C=C\\C', log)
    assert format(mol, '') == format(read_smiles(r'C/C=C\C'), '')
    assert format(mol, '') != format(read_smiles(r'C/C=C/C'), '')
    assert len(log) == 1 and 'written twice' in log[0]
    assert format(read_smiles(r'C//C=C/C'), '') == format(read_smiles(r'C/C=C/C'), '')
    with raises(IncorrectSmiles, match='two bond tokens in a row'):
        read_smiles(r'C/C=C\\\C')


def test_a_bracket_hydrogen_count_is_never_second_guessed():
    # a bracket is the writer speaking; `[CH1]` on a carbon with no bonds is wrong chemistry and is
    # stored anyway, because a reader that corrects it makes the input unrecoverable
    assert hydrogens(read_smiles('[CH1]')) == [1]


# ---------------------------------------------------------------- implicit hydrogens


def test_implicit_hydrogens_from_the_notation_model():
    cases = (
        ('C', [4]),
        ('CC', [3, 3]),
        ('CCO', [3, 2, 1]),
        ('C=C', [2, 2]),
        ('C#C', [1, 1]),
        ('N', [3]),
        ('O', [2]),
        ('S', [2]),
        ('P', [3]),
        ('F', [1]),
        ('Cl', [1]),
        ('Br', [1]),
        ('I', [1]),
        ('B', [3]),
        ('C(=O)O', [1, 0, 1]),
        # the wide valence model: N reaches 5 and I reaches 7, so a hypervalent atom gets zero
        # rather than a negative count
        ('O=N(=O)C', [0, 0, 0, 3]),
        ('S(=O)(=O)(O)O', [0, 0, 0, 1, 1]),
        ('O=I(=O)(=O)O', [0, 0, 0, 0, 1]),
    )
    for text, expected in cases:
        assert hydrogens(read_smiles(text)) == expected, text


def test_aromatic_hydrogens_come_from_the_classifier_not_from_a_second_guess():
    """Thiophene is why: charge its sulfur an extra order and S{2,4,6} lands on 4 and invents a
    hydrogen, and because that "found a valence" no fallback would ever run.  Asking
    `arom_classify_atom` first is one rule that gets all of these right.
    """
    cases = (
        ('c1ccccc1', [1] * 6),                       # benzene
        ('Cc1ccccc1', [3, 0, 1, 1, 1, 1, 1]),        # toluene: the substituted ring carbon
        ('c1ccc2ccccc2c1', [1, 1, 1, 0, 1, 1, 1, 1, 0, 1]),   # naphthalene fusion carbons
        ('c1ccncc1', [1, 1, 1, 0, 1, 1]),            # pyridine
        ('c1ccoc1', [1, 1, 1, 0, 1]),                # furan
        ('c1ccsc1', [1, 1, 1, 0, 1]),                # thiophene
        ('Cn1cccc1', [3, 0, 1, 1, 1, 1]),            # N-methylpyrrole
        ('c1cc[nH]c1', [1, 1, 1, 1, 1]),             # pyrrole: the bracket states the NH
        ('[nH0]1ccccc1', [0, 1, 1, 1, 1, 1]),        # pyridine, stated -- chython 2 cannot read it
    )
    for text, expected in cases:
        assert hydrogens(read_smiles(text)) == expected, text


def test_an_unbracketed_atom_no_rule_covers_gets_no_hydrogen_count_rather_than_zero():
    """An unbracketed atom does not state its count -- the notation's promise is that the valence
    model knows this one.  When no rule answers, 0 would be this reader inventing "none", and a
    consumer cannot tell an invented 0 from a measured one.  It can tell `None`.

    The bracketed twin is the control: a bracket DOES state the count, so its 0 is the string's own
    and must survive.  Both halves are needed -- a test that only checked the sentinel would pass
    just as well if the reader had started returning `None` for everything.
    """
    log = []
    assert hydrogens(read_smiles('N(F)(F)(F)F', log)) == [None, 0, 0, 0, 0]
    assert any('stored as unknown' in line for line in log), log

    assert hydrogens(read_smiles('[N](F)(F)(F)F')) == [0, 0, 0, 0, 0]
    assert hydrogens(read_smiles('CCO')) == [3, 2, 1]


# ---------------------------------------------------------------- what the string said


def test_aromatic_input_is_stored_aromatic():
    assert set(orders(read_smiles('c1ccccc1')).values()) == {4}


def test_the_reader_never_kekulises():
    # `changed` is False on a second call, so a True here is proof the reader did not run it first.
    # This is the invariant the whole design rests on: the caller decides when to convert.
    mol = read_smiles('c1ccccc1')
    assert kekule(mol).changed
    assert kekule(mol).changed is False


def test_kekule_input_is_stored_exactly_as_written():
    assert orders(read_smiles('C1=CC=CC=C1')) == {(1, 2): 2, (2, 3): 1, (3, 4): 2, (4, 5): 1,
                                                  (5, 6): 2, (1, 6): 1}


def test_an_explicit_aromatic_bond_between_uppercase_atoms_is_believed():
    assert orders(read_smiles('C:C')) == {(1, 2): 4}


# ---------------------------------------------------------------- atom-case promotion


BENZENE_BONDS = {(1, 2): 4, (2, 3): 4, (3, 4): 4, (4, 5): 4, (5, 6): 4, (1, 6): 4}


def test_every_all_lowercase_six_ring_reads_as_benzene():
    """The ruling, verbatim: for each smallest ring whose every atom was written lowercase, every
    bond of that ring joins the pi edge set, whatever order the string wrote.

    The aromatic set comes from atom case and never from bond order, so an explicit `-` inside an
    all-lowercase ring is a preference WITHIN the pi system.  The last two strings are exactly where
    chython 2 raises instead, and raising is an answer this reader may not give.
    """
    for text in ('c1ccccc1', 'c1cccc-c1', 'c1ccccc-1', 'c1cccc-c-1', 'c1c-c-cc-c-1'):
        assert orders(read_smiles(text)) == BENZENE_BONDS, text
        assert hydrogens(read_smiles(text)) == [1] * 6, text


def test_promotion_is_logged_as_the_repair_it_is():
    log = []
    read_smiles('c1cccc-c1', log)
    assert len(log) == 1
    assert 'lowercase' in log[0] and 'stored aromatic' in log[0]
    # and a string that needed no repair says nothing
    log = []
    read_smiles('c1ccccc1', log)
    assert log == []


def test_biphenyl_inter_ring_bond_is_not_promoted():
    """The case that makes the rule safe.  The bond between the rings lies in no smallest ring, so
    no rule here can reach it -- and `c1ccc(-c2ccccc2)cc1` is what every writer in the world emits.
    """
    for text in ('c1ccccc1-c1ccccc1', 'c1ccc(-c2ccccc2)cc1'):
        log = []
        mol = read_smiles(text, log)
        assert log == [], text
        singles = []
        for pair, order in orders(mol).items():
            if order == 1:
                singles.append(pair)
        assert len(singles) == 1, text
        # its two atoms are the ones with no hydrogen, and every other bond stayed aromatic
        assert sorted(hydrogens(mol)) == [0, 0] + [1] * 10, text


def test_biphenylene_promotes_and_still_has_a_kekule_form():
    # the four-ring joining the two six-rings is all-lowercase, so its two written-single bonds are
    # promoted; the ruling requires that the result be kekulisable, and it is
    log = []
    mol = read_smiles('c1ccc2c(c1)-c1ccccc1-2', log)
    assert len(log) == 1
    assert set(orders(mol).values()) == {4}
    assert kekule(mol).unresolved == []


def test_promotion_falls_back_to_the_stated_orders_when_it_has_no_kekule_form():
    """Free insurance: the stated set cannot do worse than itself.

    Five lowercase carbons cannot all take a ring double bond, so promoting this ring produces a
    system with no Kekule form.  The reader keeps what the string wrote and says so.
    """
    log = []
    mol = read_smiles('c1ccc-c1', log)
    assert orders(mol) == {(1, 2): 4, (2, 3): 4, (3, 4): 4, (4, 5): 1, (1, 5): 4}
    assert len(log) == 1
    assert 'no Kekule form' in log[0]


def test_promotion_redecides_every_hydrogen_count_and_rolls_them_back_with_the_orders():
    """A promoted bond changes the aromatic count of its two atoms, hence their classification,
    hence their hydrogens -- so the repair re-derives all of them rather than patching locally.

    A written TRIPLE is what makes that visible: it spends its atoms' pi electrons, so before
    promotion they are must-not-match atoms with no hydrogen, and after it they are ordinary
    aromatic CH.  Promoting a written SINGLE happens to leave carbon's bond-order sum alone, which
    is why the four benzene spellings above cannot see this at all.
    """
    # six-ring: the promotion sticks, and the two atoms of the triple gain their hydrogen
    assert hydrogens(read_smiles('c1cccc#c1')) == [1] * 6
    assert set(orders(read_smiles('c1cccc#c1')).values()) == {4}
    # five-ring: the promotion is rolled back, and so are the hydrogens.  A 1 here would mean the
    # molecule kept counts describing a graph it no longer holds.
    mol = read_smiles('c1ccc#c1')
    assert orders(mol)[(4, 5)] == 3
    assert hydrogens(mol) == [1, 1, 1, 0, 0]


def test_a_lowercase_atom_with_no_aromatic_bond_is_reported():
    log = []
    mol = read_smiles('c-c', log)
    assert orders(mol) == {(1, 2): 1}
    assert len(log) == 2
    assert 'written lowercase but carries no aromatic bond' in log[0]


# ---------------------------------------------------------------- structure


def test_ring_labels():
    assert orders(read_smiles('C1CC1')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1}
    assert orders(read_smiles('C%10CC%10')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1}
    # a label is free again once it closes
    assert orders(read_smiles('C1CC1C1CC1')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1, (3, 4): 1,
                                                 (4, 5): 1, (5, 6): 1, (4, 6): 1}


def test_a_bracketed_ring_label_reads_and_is_freed_on_close():
    # `%(NNNNN)` is the ChemAxon spelling for a label above 99. Read only: the writer numbers its own
    # closures, so nothing here comes back out as `%(...)`.
    assert orders(read_smiles('C%(101)CC%(101)')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1}
    assert orders(read_smiles('C%(0)CC%(0)')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1}
    assert orders(read_smiles('C%(99999)CC%(99999)')) == {(1, 2): 1, (2, 3): 1, (1, 3): 1}
    # a slot is per label live at once, so one string may name more of them than there are slots
    text = ''.join(f'C%({100 + i})CC%({100 + i})' for i in range(64))
    assert read_smiles(text).atom_count == 192
    # and two plain labels beside two bracketed ones do not collide
    assert orders(read_smiles('C%(101)1CC1C%(101)')) == {(1, 2): 1, (2, 3): 1, (3, 4): 1, (1, 3): 1,
                                                        (1, 4): 1}


def test_more_bracketed_labels_open_at_once_than_there_are_slots():
    with raises(IncorrectSmiles, match='open at once'):
        read_smiles(''.join(f'C%({100 + i})' for i in range(33)))


def test_ring_label_order_may_be_stated_at_either_end():
    for text in ('C=1CCCCC1', 'C1CCCCC=1'):
        assert orders(read_smiles(text))[(1, 6)] == 2, text


def test_a_ring_label_that_contradicts_itself_keeps_the_opening_order_and_says_so():
    log = []
    mol = read_smiles('C=1CCCCC-1', log)
    assert orders(mol)[(1, 6)] == 2
    assert len(log) == 1
    assert 'the opening order is kept' in log[0]


def test_branches_and_components():
    assert orders(read_smiles('CC(C)(C)C')) == {(1, 2): 1, (2, 3): 1, (2, 4): 1, (2, 5): 1}
    mol = read_smiles('CC.OO')
    assert orders(mol) == {(1, 2): 1, (3, 4): 1}
    assert len(mol.connected_components) == 2


def test_atom_numbering_follows_the_string():
    # a stable id names a token, which is what makes every log line above readable
    assert elements(read_smiles('OCN')) == [8, 6, 7]
    assert elements(read_smiles('O(C)N')) == [8, 6, 7]


# ---------------------------------------------------------------- errors


def test_syntax_errors_carry_an_offset():
    cases = (
        ('', 'no atoms'),
        ('X', "unexpected 'X' at position 0"),
        ('CCXC', "unexpected 'X' at position 2"),
        ('C(', 'unbalanced `(`'),
        ('C)', 'unbalanced `)` at position 1'),
        ('CC(C))', 'unbalanced `)` at position 5'),
        ('C1CC', 'ring bond 1 opens at position 1 and never closes'),
        ('C11', 'ring bond 1 closes on its own atom at position 2'),
        ('C12CC12', 'bonded twice'),
        ('1CC1', 'ring bond label before any atom at position 0'),
        ('C%1CC%1', '`%` needs two digits at position 1'),
        ('C%()C', '`%(` needs one to five digits and a `)` at position 1'),
        ('C%(123456)C', '`%(` needs one to five digits and a `)` at position 1'),
        ('C%(12', '`%(` needs one to five digits and a `)` at position 1'),
        ('C%(101)CC', 'ring bond 101 opens at position 1 and never closes'),
        ('[C', 'unterminated bracket atom at position 0'),
        # a bracket whose symbol position holds no letter at all: `[Zz]` is a LABEL, not an error
        ('[+2]', 'unknown element symbol at position 1'),
        ('[C:]', 'atom map `:` with no number at position 3'),
        # the top nibble value is reserved for "count unknown", so a bracket may state 14 at most
        ('[CH20]', 'hydrogen count 20 is above the storable 14'),
        ('[CH15]', 'hydrogen count 15 is above the storable 14'),
        ('[CH2H2]', 'hydrogen count given twice'),
        ('C=', 'bond token at the end of the string'),
        # the offset names the TOKEN, not where the parser noticed: `=` is at 0 and the atom that
        # made it an error is at 1
        ('=C', 'bond token starts a component at position 0'),
        ('CC.=C', 'bond token starts a component at position 3'),
        ('C=(C)C', 'bond token immediately before `(`'),
        ('C==C', 'two bond tokens in a row at position 2'),
    )
    for text, message in cases:
        with raises(IncorrectSmiles) as info:
            read_smiles(text)
        assert message in str(info.value), text


def test_tokens_this_core_refuses_by_name():
    # each of these parses somewhere else and means something this arena cannot hold; the message
    # says which, so a caller is not left guessing whether it was a typo
    # `*` is not among them: it reads as the marker, `test_smiles_r.py`
    cases = (('C$C', 'quadruple'), ('C>C', 'reaction SMILES'),
             ('C<C', 'dative bond `<-`'))
    for text, message in cases:
        with raises(IncorrectSmiles) as info:
            read_smiles(text)
        assert message in str(info.value), text


# ---------------------------------------------------------------- the dative bond


def test_the_core_can_read_its_own_dative_output():
    """`~` is order 8 here, not a SMARTS any-bond, because that is what this core's writer emits.

    Found by measurement, not by reading the spec: ammonia-borane went out of `write_smiles` as
    `[BH3]~[NH3]` and came back as an `IncorrectSmiles`.  A reader that rejects its own writer is not
    permissive-but-honest, it is broken, and no amount of "`~` means any-bond in SMARTS" fixes it.
    """
    mol = read_smiles('[BH3]~[NH3]')
    assert [b.order for b in mol.bonds()] == [8]
    assert write_smiles(mol) == '[BH3]~[NH3]'
    assert [b.order for b in read_smiles(write_smiles(mol)).bonds()] == [8]


def test_a_dative_arrow_is_read_and_the_lost_direction_is_reported():
    # `->` and `<-` are how RDKit spells a dative bond, and one turned up in a public 5k corpus
    # where it cost a whole molecule: chython 2 refuses it too, with a message about reaction SMILES
    for text in ('[NH3]->[BH3]', '[BH3]<-[NH3]'):
        log = []
        mol = read_smiles(text, log)
        assert [b.order for b in mol.bonds()] == [8], text
        # the order is storable and the arrow is not, so the arrow is reported rather than dropped
        assert any('which atom donates is not' in line for line in log), (text, log)


def test_a_dative_contact_carries_no_electron_pair():
    # ammonia donating into a metal still has three hydrogens: order 8 is out of the bond-order sum
    # and out of the aromatic classifier's neighbour count, which is what `_valence.pxi` says too
    assert hydrogens(read_smiles('N->[Fe]')) == [3, 0]
    assert hydrogens(read_smiles('N~[Fe]')) == [3, 0]
    assert hydrogens(read_smiles('C~C')) == [4, 4]


def test_chemistry_does_not_raise():
    # every one of these is impossible and every one of them comes back as a molecule with a log
    for text in ('[CH4+4]', 'C(C)(C)(C)(C)C', 'FF', '[Xe]c1ccccc1', 'c1cc[te]c1'):
        log = []
        assert read_smiles(text, log) is not None, text


# ---------------------------------------------------------------- input handling


def test_str_and_bytes_and_whitespace():
    assert elements(read_smiles('CCO')) == [6, 6, 8]
    assert elements(read_smiles(b'CCO')) == [6, 6, 8]
    assert elements(read_smiles('  CCO\n')) == [6, 6, 8]
    with raises(IncorrectSmiles):
        read_smiles('CCØ')
    with raises(TypeError):
        read_smiles(42)


def test_a_tail_is_reported_and_not_silently_dropped():
    log = []
    read_smiles('CCO ethanol', log)
    assert len(log) == 1 and 'ignored' in log[0]
    # an unterminated block is not treated as a name: it says what it is
    log = []
    read_smiles('CC |^1:0', log)
    assert len(log) == 1 and 'not terminated' in log[0]
    # the same for the brace spelling: an unclosed one is not a molecule's name
    log = []
    read_smiles('C[C@H](N)C {a:1', log)
    assert len(log) == 1 and 'not terminated' in log[0]


# ---------------------------------------------------------------- the CXSMILES tail


def radicals(mol):
    out = []
    for a in mol.atoms():
        out.append(a.is_radical)
    return out


def test_a_cx_radical_costs_its_atom_one_hydrogen():
    """`CC |^1:0|` is the ethyl radical, so it is CH2(.)-CH3 and not ethane wearing a flag.

    Measured, not asserted: the chemistry collection in `_valence.pxi` gives the same counts --
    `valence_implicit_h(z, 0, True, k) == valence_implicit_h(z, 0, False, k + 1)` for these elements
    -- and RDKit reads every one of these strings the same way.  The notation model reproducing the
    chemistry model without either file importing the other is the point.
    """
    cases = (
        ('C |^1:0|', [3]),
        ('CC |^1:0|', [2, 3]),
        ('N |^1:0|', [2]),
        ('O |^1:0|', [1]),
        # an aromatic radical too: the phenyl radical's carbon takes its ring double bond AND its
        # unpaired electron, so it has no hydrogen while the other five keep theirs
        ('c1ccccc1 |^1:0|', [0, 1, 1, 1, 1, 1]),
    )
    for text, expected in cases:
        log = []
        assert hydrogens(read_smiles(text, log)) == expected, text
        assert log == [], text                       # reading the tail is not a repair
        assert radicals(read_smiles(text))[0] is True, text


def test_a_bracket_atom_keeps_its_stated_count_when_the_tail_marks_it():
    # the bracket already answered the hydrogen question; the tail only adds the radical
    mol = read_smiles('[O]O[O] |^1:0,2|')
    assert hydrogens(mol) == [0, 0, 0]
    assert radicals(mol) == [True, False, True]
    mol = read_smiles('[CH3] |^1:0|')
    assert hydrogens(mol) == [3] and radicals(mol) == [True]


def test_a_multi_electron_radical_class_is_narrowed_and_says_so():
    # the arena has one radical bit, so a carbene cannot be stored as one.  Charging two units of
    # valence against a one-bit flag would return a molecule whose hydrogen count contradicts its own
    # radical state, which nothing downstream could tell from a real monoradical.
    log = []
    mol = read_smiles('CC |^3:0|', log)
    assert hydrogens(mol) == [2, 3] and radicals(mol) == [True, False]
    assert len(log) == 1 and 'monoradical' in log[0]


def test_the_tail_names_a_label_for_each_atom_and_they_become_aliases():
    # ChemAxon's own spelling, measured against Marvin 25.1.3: `CCC` labelled `OMe` on its last atom
    # is `CCC |$;;OMe$|`, one `;`-separated entry per atom in the tail's index space.  The label is
    # display text and the element is the file's own statement, so a labelled CARBON stays carbon --
    # unlike `[OMe]C`, where the bracket named no element and the atom is the marker.
    mol = read_smiles('CCC |$Me;;OMe$|')
    sids = [a.n for a in mol.atoms()]
    assert mol.aliases == {sids[0]: b'Me', sids[2]: b'OMe'}
    assert [a.element for a in mol.atoms()] == [6, 6, 6]
    # a label for an atom the string does not have is dropped, and the field before it still applied
    log = []
    mol = read_smiles('CC |$a;b;c$|', log)
    assert mol.aliases == {a.n: text for a, text in zip(mol.atoms(), (b'a', b'b'))}
    assert [line.rule for line in log] == ['smiles:cx-label-bad-index']


def test_a_tail_label_outranks_a_bracket_label_on_the_same_atom():
    # both name the same thing and the tail is written second, by the same producer
    log = []
    mol = read_smiles('[Pol]CC |$Resin;;$|', log)
    atom = next(a for a in mol.atoms() if a.is_r)
    assert mol.aliases == {atom.n: b'Resin'}
    assert [line.rule for line in log] == ['smiles:label-as-marker', 'smiles:cx-label-outranks-bracket']
    assert log[1].severity == 'info'


def test_a_label_spells_what_it_cannot_hold_as_a_character_reference():
    # `;` ends an entry, `$` the field, `|` the block and `&` a reference, so those four and every
    # non-ASCII character travel as `&#NN;`.  Marvin 25.1.3 writes `a;b|c$d` exactly this way.
    mol = read_smiles('CCC |$a&#59;b&#124;c&#36;d;;OMe$|')
    sids = [a.n for a in mol.atoms()]
    assert mol.aliases == {sids[0]: b'a;b|c$d', sids[2]: b'OMe'}
    # a reference's own `;` is not an entry separator, which is what makes the first assertion hold
    assert len(read_smiles('CCC |$a&#59;b;;$|').aliases) == 1
    # decimal or hexadecimal, and a code point above 127 is stored as its UTF-8 bytes
    for text in ('CC |$&#945;&#946;;$|', 'CC |$&#x3b1;&#X3B2;;$|'):
        mol = read_smiles(text)
        assert mol.aliases == {next(mol.atoms()).n: 'αβ'.encode()}, text
    # and one that is not a number, is empty or is unterminated is the text it looks like
    for text, first in (('CC |$a&#zz;b;$|', b'a&#zz'), ('CC |$&#;x;$|', b'&#'), ('CC |$a&#59;$|', b'a;')):
        mol = read_smiles(text)
        assert mol.aliases[next(mol.atoms()).n] == first, text


def test_the_reserved_labels_are_not_text():
    # `_R<n>` is ChemAxon's R-group spelling: the body writes `*` and the index is in the tail
    mol = read_smiles('CC* |$;;_R7$|')
    atom = next(a for a in mol.atoms() if a.is_r)
    assert (atom.element, atom.r_index) == (0, 7) and not mol.aliases
    # `_AP<n>` and `star_e` say what the body already says; both store nothing and say so
    for text, rule in (('CC* |$;;_AP1$|', 'smiles:cx-label-attachment-point'),
                       ('CC* |$;;star_e$|', 'smiles:cx-label-star')):
        log = []
        mol = read_smiles(text, log)
        assert not mol.aliases, text
        assert next(a for a in mol.atoms() if a.is_r).r_index == 0, text
        line = next(line for line in log if line.rule == rule)
        assert line.severity == 'info', text
    # an index on an atom the body wrote as an element has nowhere to go, and one past the domain
    # leaves the marker unindexed
    log = []
    assert read_smiles('CC |$;_R1$|', log).aliases == {}
    assert [line.rule for line in log] == ['smiles:cx-label-r-on-element']
    log = []
    mol = read_smiles('C[*] |$;_R100$|', log)
    assert next(a for a in mol.atoms() if a.is_r).r_index == 0
    assert [line.rule for line in log] == ['smiles:cx-label-r-index-too-wide']
    # both spellings at once, disagreeing: the tail is written second and wins, with a line saying so
    log = []
    mol = read_smiles('[R2]C |$_R5;$|', log)
    assert next(a for a in mol.atoms() if a.is_r).r_index == 5
    assert [line.rule for line in log] == ['smiles:cx-label-r-index-conflict']


def test_every_field_the_tail_carries_is_either_applied_or_named():
    cases = (
        ('CC |(0,0,0;1.5,0,0)|', '(...)'),                    # coordinates
        ('CC |$_AV:foo;bar$|', '$_AV:...$'),                  # atom VALUES; the labels are applied
        ('CCO |f:0.1|', 'f:0.1'),                             # fragment grouping
        ('CC |atomProp:0.p.v:1.p.w|', 'atomProp:0.p.v:1.p.w'),
    )
    for text, name in cases:
        log = []
        read_smiles(text, log)
        assert any('`%s` is not applied' % name in line for line in log), (text, log)


def test_a_field_boundary_is_found_even_after_a_field_this_reader_ignores():
    """The radical field is applied wherever it stands, which is what makes the naming safe.

    A value list may hold commas, so a comma ends a field only when a non-digit follows it; `$...$`
    and `(...)` are consumed to their closing character instead, because a `^` inside a free-text
    label is not a radical field.
    """
    for text in ('CC |f:0.1,^1:1|', 'CC |Sg:n:0,1:x:ht,^1:1|', 'CC |(0,0,0;1,0,0),^1:1|',
                 'CC |$a;b$,^1:1|', 'CC |atomProp:0.p.v,^1:1|'):
        assert radicals(read_smiles(text)) == [False, True], text
    # and a `^` inside a label block is text, not a field
    assert radicals(read_smiles('CC |$a^1:0;b$|')) == [False, False]
    # the rule cuts both ways and is stated here rather than discovered later: `,x` is indistinguish-
    # able from the start of a field called `x`, so the index before it is applied and the `x` is
    # reported as the unknown field it looks exactly like
    log = []
    assert radicals(read_smiles('CC |^1:0,x|', log)) == [True, False]
    assert any('`x` is not applied' in line for line in log), log


def test_a_broken_tail_is_reported_and_never_raises():
    # the molecule in front of the tail is intact, and refusing it would be the larger loss.  This
    # is the one place the file's "syntax raises" rule does not reach.
    cases = (
        ('CC |^1:5|', 'names atom 5, but the string has 2 atom(s)'),
        ('CC |^:0|', 'malformed'),
        ('CC |^1:|', 'names no atom'),
        ('CC |^1:0,0|', 'twice'),
        ('CC |^0:0|', 'not one of 1 to 7'),
        ('CC |^1:0x|', 'malformed after 1 index'),
    )
    for text, message in cases:
        log = []
        mol = read_smiles(text, log)
        assert mol is not None, text
        assert any(message in line for line in log), (text, log)


# ---------------------------------------------------------------- configuration
#
# WHY SOME OF THESE COMPARE RAW PARITIES AND SOME GO THROUGH INCHI.  A parity is stored against the
# unit's own `refs` order, which is built from the arena's atom order, which is the order the STRING
# created the atoms in.  So two parities are comparable only when the two strings happen to give
# their units the same frame, and each test that compares them says why it may.  Where they do not --
# any two spellings that name the atoms in a different order -- the comparison has to go through
# something canonical, and InChI is the only canonical form in reach.
#
# THE SIGNS ARE MEASURED, not adopted.  Every tetrahedral and cis/trans string below was read by
# RDKit 2026.03.4 and exported to InChI by it, and all 41 agreed with ours (2026-09-03).  Hand-written
# strings only prove the cases someone thought of, so the same comparison was then run over corpora:
# 254 configured strings carrying 1064 stated tetrahedral parities, and 1239 strings carrying a
# directional bond.  Zero opposite signs in either.  The comparison is per centre and INCLUDES the
# `/m` mirror flag, which matters more than it looks -- two enantiomers of a one-centre molecule carry
# the same `/t1-` and differ only in `/m0` against `/m1`, so a `/t`-only comparison is blind to the
# commonest case there is.  The 34 remaining mismatches all run one way, ours stating a sign where
# RDKit says undefined, on centres RDKit does not consider stereogenic; that is a perception difference
# and not this reader's to settle.  The allene could NOT be measured that way -- RDKit discards allene
# configuration on input, so `[C@]` and `[C@@]` come back as one molecule and its own SMILES output
# drops the tag -- so the allene rests on two things instead: OpenSMILES' statement that the allene
# rule IS the tetrahedral rule read over the two ends' substituents, and chython 2, which agreed with
# us on the isomer partition of twelve spellings.  The one place V2 and this reader disagree is named
# in its own test below.

def configured(mol):
    """[(anchor, parity)] over the units that carry one, so a test can state what was applied."""
    out = []
    for u in mol.stereo_units():
        if u['parity']:
            out.append((u['anchor'], u['parity']))
    return out


def frame(mol):
    """[(anchor, refs, parity)] -- the whole stored statement, for the tests that read the frame."""
    out = []
    for u in mol.stereo_units():
        out.append((u['anchor'], u['refs'], u['parity']))
    return out


def test_a_tetrahedral_configuration_is_read_over_the_written_neighbour_order():
    # `[C@H](F)(Cl)Br` and `F[C@H](Cl)Br` name their four directions in DIFFERENT orders -- the
    # bracket hydrogen counts at the position it is written, so moving it past the fluorine is one
    # transposition -- and the same tag over two orders one transposition apart is two
    # configurations.  Their frames coincide: both units order (F, Cl, Br, H), the first as
    # refs (2, 3, 4, None) and the second as (1, 3, 4, None), so the parities may be compared.
    assert frame(read_smiles('[C@H](F)(Cl)Br')) == [(1, (2, 3, 4, None), 1)]
    assert frame(read_smiles('F[C@H](Cl)Br')) == [(2, (1, 3, 4, None), 2)]
    assert frame(read_smiles('F[C@@H](Cl)Br')) == [(2, (1, 3, 4, None), 1)]


def test_the_tag_is_the_only_difference_between_two_enantiomers():
    # one spelling, two tags: the atom order is identical, so this is the one comparison that needs
    # no argument about frames at all
    for s in ('N[C@H](C)C(O)=O', 'O[C@H]1CCCC[C@H]1O', 'C[S@](=O)c1ccccc1',
              'FC(Br)=[C@]=C(F)Br', 'NC(Br)=C=[C@]=C=C(O)C'):
        a = configured(read_smiles(s))
        b = configured(read_smiles(s.replace('@', '@@', 1)))
        assert a and b and len(a) == len(b), s
        assert [p for _, p in a] != [p for _, p in b], s


@needs_inchi
def test_the_tetrahedral_sign_is_the_one_the_world_uses():
    # L-alanine's published standard InChI, and the mirror image differs only in `/m`
    l_alanine = 'InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1'
    assert molecule_to_inchi(read_smiles('C[C@@H](C(=O)O)N')) == l_alanine
    assert molecule_to_inchi(read_smiles('N[C@@H](C)C(O)=O')) == l_alanine
    assert molecule_to_inchi(read_smiles('N[C@H](C)C(O)=O')) == l_alanine.replace('/m0/', '/m1/')


@needs_inchi
def test_every_spelling_of_one_centre_is_the_same_molecule():
    # ten orders of the same four directions, and the tag chosen so each is L-alanine.  This is the
    # test that would fail if the reader read the written order approximately rather than exactly:
    # each entry differs from the one above it by at least one transposition.
    l_alanine = molecule_to_inchi(read_smiles('N[C@@H](C)C(O)=O'))
    for s in ('[C@H](N)(C)C(O)=O', 'C[C@H](N)C(O)=O', 'OC(=O)[C@@H](N)C',
              '[C@@H](C)(N)C(=O)O', 'C[C@@H](C(=O)O)N', 'N[C@H](C(O)=O)C'):
        assert molecule_to_inchi(read_smiles(s)) == l_alanine, s


@needs_inchi
def test_a_configuration_survives_a_round_trip_through_this_core_s_own_writer():
    # the anti-drift test for the sign: `translate_parity` is XOR, so the writer's parity -> tag and
    # the reader's tag -> parity are one formula read in two directions, and nothing but this notices
    # if one of them grows a correction the other does not.
    #
    # KEKULE INPUT ONLY, for now.  The writer on this branch does not write order-4 bonds -- benzene
    # comes back out as `[CH]1[CH][CH][CH][CH][CH]1`, a cyclohexane whose hydrogen counts differ --
    # so an aromatic molecule fails this for a reason that has nothing to do with configuration.
    # `C[S@](=O)c1ccccc1` is the case to add here when that lands.
    for s in ('N[C@@H](C)C(O)=O', '[C@H](F)(Cl)Br', 'F[C@H](Cl)Br', 'O[C@@H]1CCCC[C@@H]1O',
              'O[C@@H]1CCCC[C@H]1O', 'C[C@H](N)[C@@H](O)CC', 'C[S@](=O)(=O)CC',
              '[C@@H]1(O)CC[C@@H](Cl)CC1'):
        mol = read_smiles(s)
        assert molecule_to_inchi(read_smiles(write_smiles(mol))) == molecule_to_inchi(mol), s


def test_a_double_bond_reads_its_configuration_from_a_pair_of_directions():
    # `C/C=C/C` orders (C1, -, C4, -) anchored on atom 2 and `C(/C)=C/C` orders (C2, -, C4, -)
    # anchored on atom 1: the same frame shape, so the parities compare.  The branch moves the
    # direction onto the other substituent of that terminal, which is one within-pair transposition
    # and therefore inverts it -- ruling F56, and the reason a reader may not simply count slashes.
    assert frame(read_smiles('C/C=C/C'))[-1] == (2, (1, None, 4, None), 1)
    assert frame(read_smiles('C/C=C\\C'))[-1] == (2, (1, None, 4, None), 2)
    assert frame(read_smiles('C(/C)=C/C'))[-1] == (1, (2, None, 4, None), 2)
    assert frame(read_smiles('C(/C)=C\\C'))[-1] == (1, (2, None, 4, None), 1)
    # `\` on both sides is `/` on both sides upside down, which is not a different molecule
    assert frame(read_smiles('C\\C=C\\C'))[-1] == (2, (1, None, 4, None), 1)


def test_parity_one_means_the_first_and_third_directions_are_trans():
    # This reader did not choose that; `_inchi.pxi`'s `ICH_CIS_TRANS_FLIP` did, against the published
    # standard InChI of (E)-but-2-ene, and `test_inchi.py` holds the absolute assertion.  The line
    # here is the one that would notice if this reader stopped agreeing with it.
    assert configured(read_smiles('C/C=C/C')) == [(2, 1)]      # (E), methyls trans
    assert configured(read_smiles('C/C=C\\C')) == [(2, 2)]     # (Z)


@needs_inchi
def test_the_double_bond_sign_is_the_one_the_world_uses():
    assert molecule_to_inchi(read_smiles('C/C=C/C')).endswith('/b4-3+')       # (E)
    assert molecule_to_inchi(read_smiles('C/C=C\\C')).endswith('/b4-3-')      # (Z)
    # two double bonds, read independently.  A `\` shared between them is one statement about each,
    # which is what makes the middle spelling (2Z,4Z) rather than (2Z,4E)
    assert molecule_to_inchi(read_smiles('C/C=C/C=C/C')).endswith('/b5-3+,6-4+')     # (2E,4E)
    assert molecule_to_inchi(read_smiles('C/C=C\\C=C/C')).endswith('/b5-3-,6-4-')    # (2Z,4Z)
    assert molecule_to_inchi(read_smiles('C/C=C\\C=C\\C')).endswith('/b5-3-,6-4+')   # (2Z,4E)


def test_an_allene_configuration_is_read_over_the_two_ends_substituents():
    # OpenSMILES calls this allene-like and means it literally: the four substituents of the two
    # chain ends stand in for the centre's own neighbours, so the tetrahedral sentence is read over
    # them unchanged.  The centre carries the parity, and a longer odd cumulene is the same unit
    # with a longer walk to the ends.
    assert frame(read_smiles('NC(Br)=[C@]=C(O)C'))[-1] == (4, (1, 3, 6, 7), 2)
    assert frame(read_smiles('NC(Br)=[C@@]=C(O)C'))[-1] == (4, (1, 3, 6, 7), 1)
    assert frame(read_smiles('NC(Br)=C=[C@]=C=C(O)C'))[-1] == (5, (1, 3, 8, 9), 2)


@needs_inchi
def test_every_spelling_of_one_allene_is_the_same_molecule():
    # both ends written either way round and the chain entered from either end.  chython 2 partitions
    # these same eight spellings identically, which is the only cross-check available: RDKit discards
    # allene configuration outright, so it cannot arbitrate here (measured 2026-09-03).
    one = molecule_to_inchi(read_smiles('FC(Br)=[C@]=C(Cl)C'))
    other = molecule_to_inchi(read_smiles('FC(Br)=[C@@]=C(Cl)C'))
    assert one != other
    for s in ('BrC(F)=[C@@]=C(Cl)C', 'FC(Br)=[C@@]=C(C)Cl', 'ClC(C)=[C@@]=C(Br)F',
              'C(F)(Br)=[C@]=C(Cl)C'):
        assert molecule_to_inchi(read_smiles(s)) == one, s
    for s in ('ClC(C)=[C@]=C(Br)F', 'CC(Cl)=[C@@]=C(Br)F'):
        assert molecule_to_inchi(read_smiles(s)) == other, s


def test_a_hydrogen_at_an_allene_end_counts_where_the_position_rule_puts_it():
    """OpenSMILES states the "an implicit hydrogen occupies the position where it is written" rule for
    the chiral atom's own bracket; this reader applies the one rule at an allene END too, and to the
    unwritten hydrogen of a bare end as well as to a bracket's.

    So `F[CH]=` and `[CH](F)=` are two frames -- the hydrogen follows the bond to the atom before it,
    and an end that leads its component has no such bond -- and a bare end is the same frame as the
    bracket that spells its count out: `C(F)=` reads as `[CH](F)=` and `FC=` as `F[CH]=`.  Both
    equalities measured against CDK 2.12, which reads and writes the axial tag: it gives the first
    pair `CC(=[C@@]=CF)F` and the second `CC(=[C@]=CF)F`.

    A second rule for the unwritten one -- the position left over, so that `C(F)=` matched `F[CH]=`
    instead -- inverts the axial parity once per SMILES round trip, since the writer has only the one.
    """
    assert frame(read_smiles('F[CH]=[C@]=C(F)C'))[-1] == (3, (1, None, 5, 6), 2)
    assert frame(read_smiles('FC=[C@]=C(F)C'))[-1] == (3, (1, None, 5, 6), 2)
    assert frame(read_smiles('[CH](F)=[C@]=C(F)C'))[-1] == (3, (2, None, 5, 6), 1)
    assert frame(read_smiles('C(F)=[C@]=C(F)C'))[-1] == (3, (2, None, 5, 6), 1)


def test_a_configuration_this_reader_cannot_place_is_named_and_not_guessed():
    cases = (
        # `@?` is "there is a centre here and I do not know which way": a third state this arena
        # does not have, and "nobody said" is a different sentence
        ('[C@?H](F)(Cl)Br', 'states an unknown configuration'),
        # a tag on an atom that anchors no unit at all
        ('[C@](F)(Cl)(Br)(I)C', 'nothing here can hold one'),
        # a tag on a double-bond terminal: a real unit, but not a frame `@` describes
        ('C[C@H]=CC', 'is a cis/trans terminal, whose frame this reader does not build yet'),
        # two directions with no atom of their own cannot be told apart, so an order over them is
        # not an order
        ('F[C@H2]Cl', 'states 2 hydrogens and a configuration'),
        # a direction needs a partner on the far end before it says anything
        ('C/C=CC', 'has a direction on one side only'),
        # ... and a molecule that can hold no configuration at all makes every direction moot
        ('C1/C=C\\CCCC1', 'name no configuration this molecule can hold'),
        # both substituents of one terminal on the same side is a drawing that does not exist
        ('F/C(\\Cl)=C/Br', 'puts both of its substituents on the same side'),
    )
    for text, message in cases:
        log = []
        mol = read_smiles(text, log)
        assert mol is not None, text
        assert any(message in line for line in log), (text, log)


def test_a_direction_that_says_one_thing_twice_is_not_an_error():
    # `C(/F)(\Cl)=C/Br` states the same geometry on both substituents of its first terminal, which
    # is redundant and legal; only a contradiction is worth a line
    log = []
    assert configured(read_smiles('C(/F)(\\Cl)=C/Br', log)) == [(1, 2)]
    assert log == []


# ---------------------------------------------------------------- enhanced stereo groups

def test_the_tail_puts_an_atom_in_an_enhanced_stereo_group():
    assert read_smiles('F[C@H](Cl)Br |a:1|').stereo_groups() == {(STEREO_ABS, 0): [2]}
    assert read_smiles('F[C@H](Cl)Br |o1:1|').stereo_groups() == {(STEREO_OR, 1): [2]}
    assert read_smiles('F[C@H](Cl)Br |&3:1|').stereo_groups() == {(STEREO_AND, 3): [2]}
    # the field carries a list, and the configuration itself is unaffected by the grouping
    mol = read_smiles('F[C@H](Cl)Br |a:0,1,2|')
    assert mol.stereo_groups() == {(STEREO_ABS, 0): [1, 2, 3]}
    assert configured(mol) == [(2, 2)]


def test_a_group_the_tail_names_wrongly_is_reported_and_never_raises():
    cases = (
        # the range belongs to `set_stereo_group`, which is the method that validates it, so the
        # reader parses the number unbounded and reports that class's own refusal.  A number ABOVE the
        # range is renumbered instead -- the test below -- because it names a group; a stated 0 names
        # none, `o<n>` requiring the number, so it stays a refusal.
        ('F[C@H](Cl)Br |&0:1|', 'must have a group id in 1..63'),
        ('F[C@H](Cl)Br |a:9|', 'names atom 9, but the string has 4 atom(s)'),
        ('F[C@H](Cl)Br |a:|', 'names no atom'),
        ('F[C@H](Cl)Br |o:1|', '`o:1` is not applied'),           # `o` with no group number
        ('F[C@H](Cl)Br |o1:1x|', 'malformed after 1 index'),
        ('F[C@H](Cl)Br |a:1,o1:1|', 'two enhanced stereo groups'),
    )
    for text, message in cases:
        log = []
        mol = read_smiles(text, log)
        assert mol is not None, text
        assert any(message in line for line in log), (text, log)
    # 63 is the last one that fits, and it does
    assert read_smiles('F[C@H](Cl)Br |o63:1|').stereo_groups() == {(STEREO_OR, 63): [2]}


def test_a_group_id_above_the_range_is_renumbered_and_the_partition_is_what_survives():
    """A group id is a label: which atoms share a group is the statement, and the stored id is opaque
    (ruling F79), so an id the arena cannot hold is mapped to a free one of its own kind.  `union()`
    renumbers on the same grounds, and the CTfile reader repairs `MDLV30/STERAC1384` this way.

    The partition is what the assertions are about: one file id becomes one group however many atoms
    name it, an id the tail itself spends is not stolen, and the two kinds number independently.
    """
    log = []
    mol = read_smiles('C[C@H](N)[C@@H](O)[C@H](F)Cl |&1384:1,&1:3,&1384:5|', log)
    assert mol.stereo_groups() == {(STEREO_AND, 2): [2, 6], (STEREO_AND, 1): [4]}
    assert sum('renumbered to 2' in line for line in mol.log) == 1, mol.log

    assert read_smiles('F[C@H](Cl)Br.F[C@H](Cl)Br |&1:1,o64:5|').stereo_groups() == \
        {(STEREO_AND, 1): [2], (STEREO_OR, 1): [6]}, 'AND 1 does not block OR 1'

    # nothing free left: the arena's own refusal stands, and it costs that one group and no other
    spent = ','.join(f'&{i}:{i - 1}' for i in range(1, 64))
    mol = read_smiles('C' * 64 + f' |{spent},&99:63|')
    assert len(mol.stereo_groups()) == 63
    assert any('cannot be stored' in str(x) for x in mol.log), mol.log


# ---------------------------------------------------------------- the brace extension block

# Two real strings from a tool that writes the tail in braces, kept verbatim because they are the
# acceptance cases for this dialect and because each carries a field the other does not: the first
# states CIP descriptors beside an OR group, the second an AND group and no descriptors at all.
#
# Both are public structures.  The first is a steroid with a tetrahydropyranyl acetal; the second is
# a Boc-protected nitro-benzamide.  What matters for these tests is where their stereocentres are,
# which the assertions below name explicitly.
BRACE_OR = 'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CC[C@@]2(C)O[C@@H]5CCCCO5 {A19=r;A22=r;o1:19,22}'
BRACE_AND = ('Cc1c(cc(C(=O)NC[C@H]2CCN(C[C@@H]2O)C(=O)OC(C)(C)C)c3OCCCOc13)'
             '[N+](=O)[O-] {&1:9,14}')


def test_a_brace_block_states_the_same_stereo_groups_as_a_pipe_one():
    # the three group fields are spelled character for character the same in both dialects, so the
    # only thing that can differ is the scanning -- which is what this pins
    for a, b in (('F[C@H](Cl)Br {a:1}', 'F[C@H](Cl)Br |a:1|'),
                 ('F[C@H](Cl)Br {o1:1}', 'F[C@H](Cl)Br |o1:1|'),
                 ('F[C@H](Cl)Br {&3:1}', 'F[C@H](Cl)Br |&3:1|'),
                 ('F[C@H](Cl)Br {a:0,1,2}', 'F[C@H](Cl)Br |a:0,1,2|')):
        log = []
        assert read_smiles(a, log).stereo_groups() == read_smiles(b).stereo_groups(), a
        assert log == [], (a, log)


def test_a_brace_block_indexes_atoms_from_zero():
    # the negative control for every index in this section.  `1` is the stereocentre and `0` is the
    # fluorine, so a reader that counted from one would put the group on an atom that cannot hold a
    # configuration -- and would still return a molecule, silently
    assert read_smiles('F[C@H](Cl)Br {a:1}').stereo_groups() == {(STEREO_ABS, 0): [2]}


def test_the_brace_block_of_a_real_string_lands_on_its_stereocentres():
    log = []
    mol = read_smiles(BRACE_OR, log)
    # atoms 20 and 23 one-based are the `[C@@]` and the `[C@@H]`; nothing else in the string carries
    # a `@`, so this is the whole set of centres the OR group could correctly name
    assert mol.stereo_groups() == {(STEREO_OR, 1): [20, 23]}
    assert {n for n, _ in configured(mol)} == {20, 23}

    log = []
    mol = read_smiles(BRACE_AND, log)
    assert mol.stereo_groups() == {(STEREO_AND, 1): [10, 15]}
    assert {n for n, _ in configured(mol)} == {10, 15}
    # the AND string carries no `A` field, so it reads with nothing to report at all
    assert log == []


def test_a_brace_block_separates_its_fields_with_a_semicolon():
    mol = read_smiles('F[C@H](Cl)Br.F[C@H](Cl)Br {o1:1;o2:5}')
    assert mol.stereo_groups() == {(STEREO_OR, 1): [2], (STEREO_OR, 2): [6]}
    # a comma between fields is accepted as well -- the atom lists inside a field use commas, so a
    # block written by a converter between the two dialects can carry both separators
    assert read_smiles('F[C@H](Cl)Br.F[C@H](Cl)Br {o1:1,o2:5}').stereo_groups() == \
        mol.stereo_groups()


def test_a_brace_block_states_cip_descriptors():
    # `A<index>=<letter>`, which CXSMILES has no equivalent for.  Read from the string, stored in the
    # arena, and read back by stable id -- the whole path, not the reader's half of it.
    log = []
    mol = read_smiles(BRACE_OR, log)
    assert not log, log
    # 0-based in the field, 1-based as a stable id, and the same two atoms the `o1:` field named
    assert mol.atom_cips() == {20: 'r', 23: 'r'}


def test_a_brace_block_keeps_a_descriptor_s_case():
    # the reason the reader holds the letter itself and nothing upper-cases inward: lowercase r/s are
    # CIP's pseudo-asymmetric descriptors from the auxiliary rules, a different determination about a
    # different kind of centre.  A reader that normalised case would answer 'R' here and be wrong in a
    # way no round trip could see, because the writer would then be consistent with it.
    log = []
    mol = read_smiles('C[C@H](N)C {A0=R;A1=S;A2=r;A3=s}', log)
    assert not log, log
    assert mol.atom_cips() == {1: 'R', 2: 'S', 3: 'r', 4: 's'}
    assert mol.atom_cip_of(3) == 'r' and mol.atom_cip_of(3) != 'R'


def test_a_brace_block_states_a_descriptor_on_an_atom_with_no_stereo_bond():
    # storage records what the input said and does not ask whether it makes sense -- `A0` is a methyl
    # carbon here.  Asserted because the alternative (silently dropping a descriptor the reader could
    # not justify) is a repair, and a repair the caller cannot see is the one thing storage must not do.
    log = []
    mol = read_smiles('C[C@H](N)C {A0=R}', log)
    assert not log, log
    assert mol.atom_cips() == {1: 'R'}


def test_a_letter_that_is_not_an_atom_descriptor_is_refused_by_the_arena_not_the_reader():
    # 'E' is a bond descriptor, so it is a letter the reader deliberately does not judge: the domain is
    # declared once, in the arena, and its refusal becomes the log line.  This is the test that fails
    # if the reader ever grows its own copy of the accepted set.
    log = []
    mol = read_smiles('C[C@H](N)C {A1=E}', log)
    assert mol is not None
    assert len(log) == 1, log
    assert 'atom 2' in log[0] and 'cannot be stored' in log[0]
    # and the arena's own words, naming the domain it checked against
    assert "'E' is not a CIP descriptor for an atom" in log[0]
    assert not mol.atom_cips()


def test_the_q_descriptor_states_no_determination_and_is_not_a_loss():
    # `A<i>=q` is written for a centre whose descriptor the producer did not compute. Nothing to store
    # and nothing lost, so INFO -- and never a `set_atom_cip('q')` refusal, which is what a corpus of
    # 3,567 of these fields reported before.
    log = []
    mol = read_smiles('CCC1=Cc2ccc(cc2NC1=O)C(C)(C#N)Cc3ccc(cc3)C(C)C {A13=q}', log)
    assert not mol.atom_cips()
    assert len(log) == 1 and log[0].rule == 'smiles:cip-undetermined'
    assert log[0].severity == 'info'
    # beside descriptors that ARE determinations, only `q` is skipped
    log = []
    mol = read_smiles('C[C@H](N)[C@H](O)C {A1=R;A3=q}', log)
    assert mol.atom_cips() == {2: 'R'}
    assert len(log) == 1 and log[0].rule == 'smiles:cip-undetermined'


def test_the_cx_relative_flag_is_a_loss_only_where_it_stands_alone():
    # `r` carries no atom list. Beside an `&`/`o` group it restates the group -- the output is
    # byte-identical to the same string without it -- and alone it is the only statement that the
    # centres are relative, which this arena stores absolute.
    log = []
    mol = read_smiles('C[C@H](O)[C@H](N)C |&1:1,3,r|', log)
    assert format(mol, '') == format(read_smiles('C[C@H](O)[C@H](N)C |&1:1,3|'), '')
    assert len(log) == 1 and log[0].rule == 'smiles:cx-relative-flag-redundant'
    assert log[0].severity == 'info'
    log = []
    read_smiles('C[C@H](O)[C@H](N)C |r|', log)
    assert len(log) == 1 and log[0].rule == 'smiles:cx-relative-flag-unbacked'
    assert log[0].severity == 'lost'


def test_a_brace_block_field_the_reader_cannot_use_is_reported_and_never_raises():
    cases = (
        # the shape of the CIP field is this file's question and these are malformed shapes
        ('C[C@H](N)C {A1=rs}', '`A1=rs` is malformed'),
        ('C[C@H](N)C {A1=}', '`A1=` is malformed'),
        ('C[C@H](N)C {A=r}', '`A=r` is malformed'),
        ('C[C@H](N)C {A1r}', '`A1r` is malformed'),
        ('C[C@H](N)C {A1=1}', 'does not name a descriptor'),
        ('C[C@H](N)C {A9=r}', 'names atom 9, but the string has 4 atom(s)'),
        ('C[C@H](N)C {A1=R;A1=S}', 'two CIP descriptors'),
        # and the group fields report through the same path as the pipe dialect, in its words
        ('F[C@H](Cl)Br {&0:1}', 'must have a group id in 1..63'),
        ('F[C@H](Cl)Br {a:9}', 'names atom 9, but the string has 4 atom(s)'),
        ('F[C@H](Cl)Br {o:1}', '`o:1` is not applied'),
        # an unknown key is named rather than guessed at
        ('F[C@H](Cl)Br {Q1:1}', '`Q1:1` is not applied'),
        ('F[C@H](Cl)Br {f:0.1}', '`f:0.1` is not applied'),
    )
    for text, message in cases:
        log = []
        mol = read_smiles(text, log)
        assert mol is not None, text
        assert any(message in line for line in log), (text, log)


def test_a_brace_block_names_the_dialect_it_is_reporting_on():
    # the group fields are scanned by the code the pipe dialect uses, so without this the log would
    # tell a reader to look at a `|...|` tail that the string does not have
    log = []
    read_smiles('F[C@H](Cl)Br {a:9}', log)
    assert '`{...}`' in log[0] and 'CXSMILES' not in log[0], log
    log = []
    read_smiles('F[C@H](Cl)Br |a:9|', log)
    assert 'CXSMILES' in log[0] and '`{...}`' not in log[0], log


def test_the_log_is_optional_and_the_default_is_not_shared():
    # a mutable default would accumulate across calls; there is none, and no call may need one
    assert read_smiles('c1cccc-c1') is not None
    log = []
    read_smiles('c1cccc-c1', log)
    assert len(log) == 1
    read_smiles('c1cccc-c1', log)
    assert len(log) == 2
