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
"""Tests for the pass that turns a PDB record into a molecule through the residue table.

Molecules are compared as **structures** (``==`` on the container) and never as SMILES, since
SMILES-string identity is unsound here.  Every assertion about an absent log line is paired with a
fixture that does produce it, or it passes trivially."""
from pathlib import Path
from pytest import raises

from .test_pdb import _atom, _conect, _link
from ..pdb import build_molecule, mmcif, pdb


_DATA = Path(__file__).resolve().parent.parent.parent.parent / 'test'

#: The six heavy atoms of a free alanine, with the element each carries in columns 77-78.  Every
#: amino-acid row names ``OXT``, so this is the residue the table describes and not a mid-chain one.
_ALA = (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('C', ' C'), ('O', ' O'), ('OXT', ' O'))


def _deck(*lines):
    """One record out of an inline legacy deck, plus the builder's log."""
    record = pdb('\n'.join(lines))[0]
    log = []
    return build_molecule(record, log=log), log


def _residue(name, chain, seq, atoms, first=1, tag='ATOM', **kwargs):
    """``ATOM`` lines for one residue, from ``(atom_name, element)`` pairs."""
    return [_atom(first + n, atom_name, name, chain, seq, element=element, tag=tag, **kwargs)
            for n, (atom_name, element) in enumerate(atoms)]


def _orders(molecule):
    """``{(low, high): order}`` over stable ids -- the structure a bond assertion is made against."""
    return {(bond.n, bond.m) if bond.n <= bond.m else (bond.m, bond.n): bond.order
            for bond in molecule.bonds()}


def _names(molecule):
    """The element of every atom, in container order, which is the file's order."""
    return [atom.atomic_symbol for atom in molecule.atoms()]


def _hydrogens(molecule):
    return [molecule.implicit_h_of(atom.n) for atom in molecule.atoms()]


def _lines(log, probe):
    return [str(x) for x in log if probe in x]


# the template, applied


def test_one_complete_alanine_gets_the_templates_five_bonds_and_its_hydrogens():
    """A free ALA whose six heavy atoms the file states comes back as alanine.

    The legacy deck states no bond, so all five are the table's -- including the ``C=O`` order 2, which
    a file of this format cannot say.  Hydrogen counts follow from those orders.  [mutant: skipping the
    `template.bonds` loop in `build_molecule`]
    """
    molecule, log = _deck(*_residue('ALA', 'A', 1, _ALA))
    assert _names(molecule) == ['N', 'C', 'C', 'C', 'O', 'O']
    # N-CA, CA-CB, CA-C, C-O (double), C-OXT
    assert _orders(molecule) == {(1, 2): 1, (2, 3): 1, (2, 4): 1, (4, 5): 2, (4, 6): 1}
    assert _hydrogens(molecule) == [2, 1, 3, 0, 0, 1]
    # armed negative for the shortfall line below: a residue with every atom its row names is silent
    assert log == []


def test_a_residue_missing_an_atom_loses_that_atoms_bonds_and_nothing_else():
    """ALA with no ``CB``: four bonds, no fifth atom invented to carry the fifth.

    Both endpoints or no bond: a residue the file drew incompletely acquires neither the bond nor the
    atom.  [mutant: dropping the `if a is None or b is None: continue` guard in `build_molecule`'s
    template-bond loop]
    """
    molecule, log = _deck(*_residue('ALA', 'A', 1, [p for p in _ALA if p[0] != 'CB']))
    assert _names(molecule) == ['N', 'C', 'C', 'O', 'O']
    assert _orders(molecule) == {(1, 2): 1, (2, 3): 1, (3, 4): 2, (3, 5): 1}

    shortfall = _lines(log, 'lack the template atom')
    assert len(shortfall) == 1, log
    assert shortfall[0] == ('residue: 1 residue(s) of ALA lack the template atom(s) CB; no bond to an '
                            'absent atom is applied and no atom is invented (first A/ALA1)')


def test_a_mid_chain_residue_missing_only_its_link_displacement_is_silent():
    """A mid-chain GLY missing only OXT produces no shortfall line.

    Every amino acid row names OXT, but a mid-chain residue legitimately lacks it: the chain link
    displaced it.  The suppression fires because OXT has exactly one template neighbour (C) and the
    link at C was built.  Armed by its twin below.  [mutant: removing the suppression loop in
    ``_shortfall``]
    """
    # GLY seq 1 (mid-chain, no OXT): N CA C O; GLY seq 2 (terminal, has OXT): N CA C O OXT
    molecule, log = _deck(*_gly(1, 1), *_residue('GLY', 'A', 2,
                                                 (('N', ' N'), ('CA', ' C'), ('C', ' C'),
                                                  ('O', ' O'), ('OXT', ' O')), first=5))
    assert not _lines(log, 'lack the template atom')
    assert molecule.bond_count == 8              # 3 + 3 + 1 peptide + 1 C-OXT in seq 2


def test_a_residue_genuinely_missing_an_atom_still_reports_it():
    """The positive twin: a mid-chain GLY that is also missing its backbone O still produces the line.

    GLY missing both O and OXT has two candidates at the link atom C; the budget of one suppresses
    neither, so both are reported.  [mutant: suppressing all candidates regardless of count in
    ``_shortfall``]
    """
    # GLY seq 1 missing both O and OXT; GLY seq 2 has all atoms including OXT (terminal, no link at C)
    gly1_no_o_oxt = _residue('GLY', 'A', 1, (('N', ' N'), ('CA', ' C'), ('C', ' C')), first=1)
    gly2_full = _residue('GLY', 'A', 2, (('N', ' N'), ('CA', ' C'), ('C', ' C'),
                                         ('O', ' O'), ('OXT', ' O')), first=5)
    molecule, log = _deck(*gly1_no_o_oxt, *gly2_full)
    shortfall = _lines(log, 'lack the template atom')
    assert len(shortfall) == 1, log
    assert shortfall[0] == ('residue: 1 residue(s) of GLY lack the template atom(s) O, OXT; no bond '
                            'to an absent atom is applied and no atom is invented (first A/GLY1)')


def test_a_mid_chain_residue_missing_two_candidates_reports_both():
    """A mid-chain ALA missing both its carbonyl O and terminal OXT reports both, not just OXT.

    A built C link accounts for exactly one absent neighbour, so two candidates means one is
    unaccounted for and both are reported.  [mutant: suppressing both candidates when there are two]
    """
    # ALA seq 1 missing O and OXT (both bonded only to C); ALA seq 2 intact for the link
    ala_no_o_oxt = _residue('ALA', 'A', 1, (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('C', ' C')),
                            first=1)
    ala2 = _residue('ALA', 'A', 2, (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('C', ' C'),
                                    ('O', ' O'), ('OXT', ' O')), first=10)
    _, log = _deck(*ala_no_o_oxt, *ala2)
    shortfall = _lines(log, 'lack the template atom')
    assert len(shortfall) == 1, log
    assert 'O, OXT' in shortfall[0], shortfall


def test_three_residues_short_of_one_atom_are_one_finding_with_a_count_of_three():
    """A protein has thousands of residues, so the shortfall is counted and not narrated.

    The count is asserted, not the mere presence of a line.  [mutant: `self.counts[key] = 1` in
    `_Findings.hit`]
    """
    lines = []
    for seq in (1, 2, 3):
        lines += _residue('ALA', 'A', seq, [p for p in _ALA if p[0] != 'CB'], first=10 * seq)
    _, log = _deck(*lines)
    shortfall = _lines(log, 'lack the template atom')
    assert len(shortfall) == 1, log
    assert shortfall[0].startswith('residue: 3 residue(s) of ALA lack the template atom(s) CB;')
    assert shortfall[0].endswith('(first A/ALA1)')


def test_an_aromatic_residue_arrives_kekule_and_thiele_aromatizes_it():
    """HIS comes off its row with alternating ring orders, and ``thiele()`` closes the loop.

    ``thiele()`` refuses rather than repairing, so its accepting this ring is the assertion: the orders
    it was handed were a valid Kekule form.  [mutant: flattening every template order to 1 in
    `build_molecule`'s template-bond loop]
    """
    his = (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('CG', ' C'), ('ND1', ' N'), ('CD2', ' C'),
           ('CE1', ' C'), ('NE2', ' N'), ('C', ' C'), ('O', ' O'), ('OXT', ' O'))
    molecule, log = _deck(*_residue('HIS', 'A', 1, his))
    assert log == []
    assert molecule.is_kekule
    assert molecule.aromatic_rings_count == 0
    assert sorted(_orders(molecule).values()) == [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]

    assert molecule.thiele()
    assert molecule.aromatic_rings_count == 1
    assert molecule.aromatic_bond_count == 5


def test_water_is_one_atom_with_two_hydrogens_and_an_ion_carries_the_tables_charge():
    """Two single-atom rows, and the charge on one of them comes from the table.

    The sodium's charge column is blank, and a blank column is not a claim of neutrality -- the record
    cannot tell it from a stated zero, so it yields to the table.  [mutant: `charge = atom.charge`
    unconditionally in `_atom_spec`]
    """
    molecule, log = _deck(_atom(1, 'O', 'HOH', 'W', 1, element=' O', tag='HETATM'),
                          _atom(2, 'NA', ' NA', 'B', 1, element='NA', tag='HETATM'))
    assert log == []
    assert _names(molecule) == ['O', 'Na']
    assert _hydrogens(molecule) == [2, 0]
    assert [atom.charge for atom in molecule.atoms()] == [0, 1]
    assert not molecule.bond_count


def test_a_non_zero_charge_from_the_file_beats_the_table_and_says_so():
    """The other half of the charge rule, and the armed negative for the line above.

    A stated non-zero charge is a claim about *this* atom and wins over one about a component.
    [mutant: `charge = row_charge` unconditionally in `_atom_spec`]
    """
    molecule, log = _deck(_atom(1, 'NA', ' NA', 'B', 1, element='NA', charge='2+', tag='HETATM'))
    assert [atom.charge for atom in molecule.atoms()] == [2]
    assert _lines(log, 'state charge 2 where the template has 1') == [
        "atom: 1 atom(s) state charge 2 where the template has 1; the file's charge is used "
        '(first B/NA1 NA)']


def test_the_files_element_beats_the_table_and_the_mismatch_is_reported():
    """A misnamed atom keeps the element the file gave it, and gets a line.

    The file's element is a statement about this atom; ours is one about a component.  [mutant:
    `element = row_element` unconditionally in `_atom_spec`]
    """
    molecule, log = _deck(*_residue('ALA', 'A', 1,
                                    [(name, ' N' if name == 'CB' else element)
                                     for name, element in _ALA]))
    assert _names(molecule) == ['N', 'C', 'N', 'C', 'O', 'O']
    assert _lines(log, 'state element N where the template has C') == [
        "atom: 1 atom(s) state element N where the template has C; the file's element is used "
        '(first A/ALA1 CB)']


# the polymer link


def _gly(seq, first, chain='A'):
    return _residue('GLY', chain, seq, (('N', ' N'), ('CA', ' C'), ('C', ' C'), ('O', ' O')),
                    first=first)


def test_two_consecutive_residues_get_exactly_one_peptide_bond():
    """``C`` of residue 1 to ``N`` of residue 2, once, and nothing else joins the two.

    The count matters as much as the bond: a rule firing per atom pair rather than per residue pair
    joins more than the two link atoms.  [mutant: dropping the `pairs.setdefault` call in
    `_chain_links`]
    """
    molecule, log = _deck(*_gly(1, 1), *_gly(2, 5))
    assert molecule.connected_components_count == 1
    # 3 bonds within each GLY (N-CA, CA-C, C=O) plus the one peptide bond
    assert molecule.bond_count == 7
    assert _orders(molecule)[(3, 5)] == 1        # GLY1 C (stable id 3) to GLY2 N (stable id 5)
    assert not _lines(log, 'chain link')


def test_a_numbering_gap_leaves_the_chain_broken_and_names_both_residues():
    """Residues 1 and 3 of one chain are two molecules, and the log says which pair.

    It cannot be checked further: a distance test would be perception and this package has none, so
    what is reported is the numbering and nothing more.  [mutant: accepting any delta in
    `_chain_links`]
    """
    molecule, log = _deck(*_gly(1, 1), *_gly(3, 5))
    assert molecule.connected_components_count == 2
    assert molecule.bond_count == 6

    gap = _lines(log, 'sequence numbers skip')
    assert len(gap) == 1, log
    assert gap[0] == ('residue: 1 residue pair(s) are adjacent in their chain but their sequence '
                      'numbers skip, so no chain link is built between them (first A/GLY1-A/GLY3)')


def test_two_chains_are_never_linked_across():
    """Residue 1 of chain A and residue 2 of chain B are two molecules and no gap.

    The chain is part of the grouping key, so the two are never adjacent in one sort -- which is why
    this produces no line either.  [mutant: dropping `key[1]` from the chain key in `_chain_links`]
    """
    molecule, log = _deck(*_gly(1, 1, 'A'), *_gly(2, 5, 'B'))
    assert molecule.connected_components_count == 2
    assert not _lines(log, 'sequence numbers skip')


def test_the_link_is_not_built_when_a_link_atom_is_absent():
    """A residue whose ``C`` the file never stated cannot be joined to the next one.

    Both endpoints or no bond, applied to the link.  Three components and not two, because ``O`` bonds
    to nothing but the missing ``C`` and falls off with it.  [mutant: removing the `if a is None or b is
    None` guard in `_chain_links`]
    """
    first = _residue('GLY', 'A', 1, (('N', ' N'), ('CA', ' C'), ('O', ' O')), first=1)
    molecule, log = _deck(*first, *_gly(2, 5))
    assert molecule.connected_components_count == 3
    assert _lines(log, 'chain link(s) are not built') == [
        'residue: 1 chain link(s) are not built because a link atom the template names (C or N) is '
        'absent (first A/GLY1-A/GLY2)']


def test_the_chain_is_sorted_and_not_taken_in_file_order():
    """The same two residues written in reverse order give the same one peptide bond.

    File order is only conventionally sequence order.  [mutant: dropping the `group.sort` call in
    `_chain_links`]
    """
    forward, _ = _deck(*_gly(1, 1), *_gly(2, 5))
    reverse, log = _deck(*_gly(2, 5), *_gly(1, 1))
    assert reverse.bond_count == 7
    assert forward == reverse
    assert not _lines(log, 'sequence numbers skip')


def test_a_stated_link_and_the_numbering_agree_on_one_bond():
    """A ``LINK`` record naming the same two atoms does not produce a second peptide bond.

    The two routes converge because both spell the pair low-first; ``LINK`` states no order, so the key
    normalisation and not the order ``setdefault`` is what the count rests on.  [mutant: `key = (b, a)`
    in `_chain_links`]
    """
    molecule, log = _deck(*_gly(1, 1), *_gly(2, 5),
                          _link('C', 'GLY', 'A', 1, 'N', 'GLY', 'A', 2))
    assert molecule.bond_count == 7
    assert _orders(molecule)[(3, 5)] == 1
    assert not _lines(log, 'sequence numbers skip')


# bond orders


def _cif(comp, atoms, rows, tag='ATOM'):
    """A minimal mmCIF stating one residue and, in ``_chem_comp_bond``, the rows given.

    ``_chem_comp_bond`` is the only route to a bond whose order the file states -- a legacy ``CONECT``
    states none -- so it is how the file-versus-template order rule is reachable at all.
    """
    atoms = '\n'.join(
        f'{tag} {n} {element.strip()} {name} . {comp} A 1 1 ? {n}.000 0.000 0.000 1.00 20.00 ? 1 A 1'
        for n, (name, element) in enumerate(atoms, 1))
    bonds = '\n'.join(f'{comp} {a} {b} {value}' for a, b, value in rows)
    return (f'data_{comp}\n'
            'loop_\n_entity.id\n_entity.type\n1 polymer\n#\n'
            'loop_\n'
            '_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n'
            '_chem_comp_bond.value_order\n'
            f'{bonds}\n#\n'
            'loop_\n'
            '_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n'
            '_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n'
            '_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n'
            '_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n'
            '_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n'
            '_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n'
            '_atom_site.pdbx_PDB_model_num\n'
            f'{atoms}\n#\n')


def _chem_comp_cif(*rows):
    """The same, for the one alanine every order test below is written against."""
    return _cif('ALA', _ALA, rows)


def _from_cif(text):
    log = []
    return build_molecule(mmcif(text)[0], log=log), log


def test_a_stated_bond_duplicating_a_template_bond_is_one_bond():
    """``_chem_comp_bond`` states the bonds the table states, and the pair is deduplicated.

    A distributed mmCIF entry states every intra-residue bond, so without deduplication by the
    unordered pair every residue comes back with each of its bonds twice.  [mutant: `pairs[key] =
    bond.order` unconditionally in `_stated_bonds`]
    """
    molecule, log = _from_cif(_chem_comp_cif(('C', 'O', 'DOUB'), ('CA', 'CB', 'SING')))
    assert molecule.bond_count == 5
    assert _orders(molecule) == {(1, 2): 1, (2, 3): 1, (2, 4): 1, (4, 5): 2, (4, 6): 1}
    assert not _lines(log, 'where the template has')


def test_a_stated_order_disagreeing_with_the_template_wins_and_says_so():
    """``SING`` on the carbonyl is the file contradicting our table, and the file is more specific.

    It made an explicit claim about *this* bond; the table makes one about a component.  [mutant:
    dropping the `pairs[key] = bond.order` assignment in the disagreement branch of `_stated_bonds`]
    """
    molecule, log = _from_cif(_chem_comp_cif(('C', 'O', 'SING')))
    assert _orders(molecule)[(4, 5)] == 1
    disagreement = _lines(log, 'where the template has')
    assert len(disagreement) == 1, log
    assert disagreement[0] == ("bond: 1 stated bond(s) state order 1 where the template has 2; the "
                               "file's order is used (first atoms 3-4)")


def test_an_unstated_order_over_a_template_double_bond_stays_double_and_says_nothing():
    """A ``CONECT`` between the carbonyl C and its O is order 2, silently.

    ``CONECT`` states no order in any writer, so the record's order 1 is convention and not evidence,
    and the template is the only evidence there is.  The silence is deliberate: it is not a
    disagreement.  [mutant: ignoring `bond.stated_order` in `_stated_bonds`]
    """
    molecule, log = _deck(*_residue('ALA', 'A', 1, _ALA), _conect(4, 5))
    assert _orders(molecule) == {(1, 2): 1, (2, 3): 1, (2, 4): 1, (4, 5): 2, (4, 6): 1}
    assert not _lines(log, 'where the template has')
    assert log == []


# hydrogens


def test_explicit_hydrogens_on_a_templated_residue_are_dropped_and_the_count_matches():
    """A file that carries alanine's seven hydrogens gives the same molecule as one that does not.

    ``ResidueTemplate.atoms`` is heavy atoms only, so an explicit hydrogen has no template bond and
    would otherwise arrive as an isolated atom.  [mutant: removing the `atom.element == 'H'` branch in
    `_plan_templated`]
    """
    protons = ('H', 'H2', 'HA', 'HB1', 'HB2', 'HB3', 'HXT')
    lines = _residue('ALA', 'A', 1, _ALA)
    lines += _residue('ALA', 'A', 1, [(name, ' H') for name in protons], first=20)
    molecule, log = _deck(*lines)
    bare, _ = _deck(*_residue('ALA', 'A', 1, _ALA))

    assert molecule.atom_count == 6
    assert molecule == bare
    assert sum(_hydrogens(molecule)) == len(protons)
    assert log == []


def test_a_protonation_the_derivation_disagrees_with_is_reported_with_a_count():
    """Six explicit hydrogens on an alanine that derives seven: the drop becomes a report.

    The armed half of the test above.  Aggregated by residue name and by the two counts, because a
    systematic protonation difference is one finding repeated.  [mutant: removing the `derived !=
    residue.stated_h` check in `_hydrogens`]
    """
    protons = ('H', 'HA', 'HB1', 'HB2', 'HB3', 'HXT')
    lines = []
    # two chains rather than two sequence numbers: consecutive numbering would build the peptide bond,
    # and a linked residue derives a different count for an unrelated reason
    for chain in ('A', 'B'):
        lines += _residue('ALA', chain, 1, _ALA, first=20 * (chain == 'B') + 1)
        lines += _residue('ALA', chain, 1, [(n, ' H') for n in protons],
                          first=20 * (chain == 'B') + 11)
    _, log = _deck(*lines)

    mismatch = _lines(log, 'explicit hydrogen(s) where')
    assert len(mismatch) == 1, log
    assert mismatch[0] == ('atom: 2 residue(s) of ALA state 6 explicit hydrogen(s) where 7 are '
                           'derived; the explicit hydrogens are dropped and the derived count is used '
                           '(first A/ALA1)')


def test_a_residue_with_an_unsettled_atom_is_not_compared_against_its_stated_h_count():
    """A residue holding an ``H_UNKNOWN`` atom and a stated explicit hydrogen count produces no mismatch.

    ``H_UNKNOWN`` means the derivation cannot answer, and an unsettled count is not a statement -- an
    ``or 0`` would understate the derived total and fire the mismatch on a fine residue.  ``AROM`` in
    ``_chem_comp_bond`` overrides the template's Kekule orders, leaving the two ring nitrogens with no
    derivable count.  Armed by the twin below.  [mutant: removing the `any(... is None ...)` skip in
    `_hydrogens`]
    """
    # HIS with its imidazole ring stated AROM, so the two ring nitrogens become H_UNKNOWN
    his_atoms = (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('CG', ' C'), ('ND1', ' N'),
                 ('CD2', ' C'), ('CE1', ' C'), ('NE2', ' N'), ('C', ' C'), ('O', ' O'),
                 ('OXT', ' O'))
    arom_bonds = [('ND1', 'CG', 'AROM'), ('CD2', 'CG', 'AROM'), ('CE1', 'ND1', 'AROM'),
                  ('NE2', 'CD2', 'AROM'), ('NE2', 'CE1', 'AROM')]
    from ..pdb import mmcif

    def _his_cif_arom():
        atom_lines = '\n'.join(
            f'ATOM {n} {el.strip()} {name} . HIS A 1 1 ? {n}.000 0.000 0.000 1.00 20.00 ? 1 A 1'
            for n, (name, el) in enumerate(his_atoms, 1))
        # one explicit H so that stated_h > 0; fields must align: type_symbol H, atom_name H
        atom_lines += (f'\nATOM {len(his_atoms)+1} H H . HIS A 1 1 ? '
                       f'{len(his_atoms)+1}.000 0.000 0.000 1.00 20.00 ? 1 A 1')
        bond_lines = '\n'.join(f'HIS {a} {b} {v}' for a, b, v in arom_bonds)
        return (f'data_HIS\nloop_\n_entity.id\n_entity.type\n1 polymer\n#\n'
                'loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n'
                '_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n'
                f'{bond_lines}\n#\n'
                'loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n'
                '_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n'
                '_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n'
                '_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n'
                '_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n'
                '_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n'
                '_atom_site.pdbx_PDB_model_num\n'
                f'{atom_lines}\n#\n')

    _, log = _from_cif(_his_cif_arom())
    assert not _lines(log, 'explicit hydrogen(s) where')


def test_a_residue_with_all_settled_atoms_and_a_wrong_stated_count_still_reports_the_mismatch():
    """The armed twin: all counts settled, wrong stated H, the mismatch line fires.

    HIS with its Kekule template orders has every ring atom settled after ``derive_hydrogens``, so one
    stated explicit H against a derived nine is a genuine discrepancy.  [mutant: removing the `derived
    != residue.stated_h` check in `_hydrogens`]
    """
    his = (('N', ' N'), ('CA', ' C'), ('CB', ' C'), ('CG', ' C'), ('ND1', ' N'), ('CD2', ' C'),
           ('CE1', ' C'), ('NE2', ' N'), ('C', ' C'), ('O', ' O'), ('OXT', ' O'))
    lines = _residue('HIS', 'A', 1, his)
    # One explicit H stated; the derivation gives nine.
    lines += _residue('HIS', 'A', 1, [('H', ' H')], first=20)
    _, log = _deck(*lines)
    mismatch = _lines(log, 'explicit hydrogen(s) where')
    assert len(mismatch) == 1, log
    assert '1 explicit hydrogen(s) where 9 are derived' in mismatch[0]


def test_a_heavy_atom_only_residue_is_not_checked_against_zero():
    """A file that states no hydrogen at all has made no statement about protonation.

    A check against zero here fires on every residue of almost every archive entry.  Armed by the test
    above.  [mutant: removing the `if not residue.stated_h` guard in `_hydrogens`]
    """
    _, log = _deck(*_residue('ALA', 'A', 1, _ALA), *_gly(2, 20))
    assert not _lines(log, 'explicit hydrogen(s) where')


# the ligand


def test_a_residue_with_no_template_keeps_its_atoms_and_only_the_stated_bonds():
    """A ligand the table does not know: every atom, one stated bond, and nothing perceived.

    The ``unsupported: `` prefix claims the file is fine and our 57-row table is the limitation; the
    unbonded count names ``saturate()`` as the next step without calling it.  The order assertion is
    the substantive half: nothing acquired an order from a search.  [mutant: dropping the
    `unsupported: ` log line from `build_molecule`]
    """
    ligand = (('C1', ' C'), ('C2', ' C'), ('O1', ' O'), ('CL', 'CL'))
    molecule, log = _deck(*_residue('LIG', 'A', 1, ligand, tag='HETATM'), _conect(1, 2))
    assert _names(molecule) == ['C', 'C', 'O', 'Cl']
    assert _orders(molecule) == {(1, 2): 1}

    unsupported = [str(x) for x in log if str(x).startswith('unsupported: ')]
    assert unsupported == ['unsupported: 1 residue(s) of 1 component id(s) have no row in the residue '
                           'table, so no template bond is applied to them: LIG (1)']
    assert _lines(log, 'hold no bond at all') == [
        'residue: 2 atom(s) of residue(s) with no template row hold no bond at all; nothing is '
        'invented for them, and chython.chemistry.saturate() is the separate pass a caller runs on a '
        'ligand whose file gave connectivity and no orders']


def test_a_templated_residue_never_reports_an_unbonded_count():
    """The armed negative for the line above: a water is one unbonded atom and is not a ligand.

    A single-atom residue the table *does* know has nothing to bond to, which is not the ``saturate()``
    population.  [mutant: counting every atom rather than the members of `loose` in `build_molecule`]
    """
    _, log = _deck(_atom(1, 'O', 'HOH', 'W', 1, element=' O', tag='HETATM'))
    assert not _lines(log, 'hold no bond at all')
    assert log == []


def test_an_unsettled_hydrogen_count_is_reported_and_kekule_is_named_as_the_repair():
    """An untemplated ligand whose ring the file states as aromatic leaves one atom unsettled.

    The one case the shared derivation does not answer: the pnictogen whose class the ring decides.
    The aromatic carbons are settled, so the count is 1 and not 5, and the line names the repair rather
    than performing it.  [mutant: `unsettled = {}` in `_hydrogens`]
    """
    ring = (('N1', ' N'), ('C2', ' C'), ('C3', ' C'), ('C4', ' C'), ('C5', ' C'))
    pairs = (('N1', 'C2'), ('C2', 'C3'), ('C3', 'C4'), ('C4', 'C5'), ('C5', 'N1'))
    molecule, log = _from_cif(_cif('LIG', ring, [(a, b, 'AROM') for a, b in pairs], tag='HETATM'))
    assert molecule.aromatic_rings_count == 1
    assert _lines(log, 'no derivable implicit hydrogen count') == [
        'atom: 1 atom(s) hold no derivable implicit hydrogen count; kekule() settles the aromatic '
        'pnictogen and check_valence() names the rest']

    # the line names one call, and the call settles it -- a repair the caller runs, never the builder
    assert molecule.kekule()
    assert not molecule.derive_hydrogens()


def test_an_atom_the_template_does_not_name_is_kept_and_reported():
    """A stated atom is never dropped for disagreeing with our table; only a hydrogen is.

    It gets no template bond, which is all the table can say about it.  [mutant: `continue` instead of
    storing the atom in the `row is None` branch of `_plan_templated`]
    """
    molecule, log = _deck(*_residue('ALA', 'A', 1, _ALA + (('SE', 'SE'),)))
    assert _names(molecule) == ['N', 'C', 'C', 'C', 'O', 'O', 'Se']
    assert molecule.bond_count == 5
    assert _lines(log, 'the ALA template does not have') == [
        'residue: 1 atom(s) of ALA carry a name the ALA template does not have (SE); they are stored '
        'and get no template bond (first A/ALA1)']


# alternate conformers


def _ser_conformers():
    """One SER whose ``CB``/``OG`` are drawn twice, conformer ``B`` at the higher occupancy."""
    shared = _residue('SER', 'A', 1, (('N', ' N'), ('CA', ' C'), ('C', ' C'), ('O', ' O')))
    side = (('CB', ' C'), ('OG', ' O'))
    # the two conformers differ only in x, which is the only way a test can tell which was taken:
    # they describe the same two atoms and give the same graph either way
    return (shared + _residue('SER', 'A', 1, side, first=5, alt='A', occupancy=0.4, x=1.)
            + _residue('SER', 'A', 1, side, first=7, alt='B', occupancy=0.6, x=2.))


def test_one_conformer_is_selected_by_occupancy_and_the_selection_is_logged():
    """Two descriptions of one atom cannot both be in the graph, so one is chosen and named.

    By summed occupancy and not by file order.  ``B`` is written second here, so "the first one seen"
    would pick ``A``.  [mutant: `max` in place of `min`, or dropping the occupancy term, in
    `_select_conformers`]
    """
    molecule, log = _deck(*_ser_conformers())
    assert molecule.atom_count == 6
    assert molecule.connected_components_count == 1

    selection = _lines(log, 'alternate conformers')
    assert len(selection) == 1, log
    assert selection[0] == ("residue: 1 residue(s) hold alternate conformers; the one with the "
                            "highest summed occupancy is kept and the rest are dropped "
                            "(first A/SER1 keeps 'B' of 'A', 'B')")


def test_file_order_does_not_decide_which_conformer_wins():
    """The same two conformers written the other way round select the same one.

    The structures are equal either way -- the graph cannot tell two conformers apart -- so the
    assertion has to be the coordinate: ``x=2.0`` is conformer ``B``.  [mutant: `min(alternates,
    key=lambda alt: by_alt[alt][0])` in `_select_conformers`]
    """
    lines = _ser_conformers()
    swapped = lines[:4] + lines[6:] + lines[4:6]
    first, _ = _deck(*lines)
    second, _ = _deck(*swapped)
    assert first == second
    assert first.xy_of(5)[0] == 2.
    assert second.xy_of(5)[0] == 2.


def test_an_explicit_alt_loc_selects_that_conformer():
    """``alt_loc='A'`` takes the lower-occupancy conformer, because the caller asked for it.

    [mutant: ignoring the `requested` argument in `_select_conformers`]
    """
    default, _ = _deck(*_ser_conformers())
    record = pdb('\n'.join(_ser_conformers()))[0]
    explicit_log = []
    explicit = build_molecule(record, alt_loc='A', log=explicit_log)
    assert explicit.atom_count == 6
    # the graph is the same either way, so the coordinate says which was taken: x=1.0 'A', x=2.0 'B'
    assert explicit.xy_of(5)[0] == 1.
    assert default.xy_of(5)[0] == 2.
    assert not _lines(explicit_log, 'alternate conformers')


def test_a_residue_lacking_the_requested_conformer_is_logged():
    """``alt_loc='C'`` on a residue that has ``A`` and ``B`` keeps only its conformer-free atoms.

    Falling back to another id would answer a question the caller did not ask.  [mutant: falling back
    to the occupancy selection when the requested id is absent]
    """
    record = pdb('\n'.join(_ser_conformers()))[0]
    log = []
    molecule = build_molecule(record, alt_loc='C', log=log)
    assert molecule.atom_count == 4
    assert _lines(log, 'not the requested') == [
        "residue: 1 residue(s) hold alternate conformers but not the requested 'C', so only their "
        "conformer-free atoms are kept (first A/SER1 holds 'A', 'B')"]


def test_a_bond_is_never_built_across_two_conformers():
    """A ``CONECT`` between conformer ``A``'s ``CB`` and conformer ``B``'s ``OG`` is not a bond.

    Two atoms with different alt ids describe the same place (see ``PDBAtom.residue_key``), and the
    selection always drops one endpoint.  The ``== 5`` assertion measures the guard only when paired
    with the test below.  [mutant: dropping the `bond.a not in kept` guard in `_stated_bonds`]
    """
    molecule, log = _deck(*_ser_conformers(), _conect(5, 8))
    assert molecule.atom_count == 6
    assert molecule.bond_count == 5              # N-CA, CA-C, C=O, CA-CB, CB-OG
    assert not _lines(log, 'could not be stored')


def test_a_bond_is_built_when_both_endpoints_survive_conformer_selection():
    """A ``CONECT`` naming two atoms of the selected conformer builds the bond.

    Armed twin of the test above.  Serial 1 is N (always kept), serial 8 is OG-B (the selected
    conformer); N-OG is not a template bond, so its only source is the CONECT.  [mutant: inverting the
    guard in `_stated_bonds` to `if bond.a in kept or bond.b in kept`]
    """
    molecule, log = _deck(*_ser_conformers(), _conect(1, 8))
    assert molecule.atom_count == 6
    assert molecule.bond_count == 6              # template 5 + N-OG from CONECT
    assert not _lines(log, 'could not be stored')


def test_a_stated_bond_naming_an_unstored_atom_is_reported_and_not_built():
    """A ``CONECT`` pointing at an atom that could not be stored produces a bond finding.

    An atom of an untemplated residue that states no element cannot be stored, so a bond to it has a
    missing endpoint.  Positive twin of ``test_a_bond_is_never_built_across_two_conformers``.  [mutant:
    removing the `bond.a not in plan or bond.b not in plan` guard in `_stated_bonds`]
    """
    # atom 1 has an element and is stored; atom 2 has none and is not, so the CONECT loses an endpoint
    molecule, log = _deck(
        _atom(1, 'C1', 'LIG', 'A', 1, element=' C', tag='HETATM'),
        _atom(2, 'X1', 'LIG', 'A', 1, element='', tag='HETATM'),
        _conect(1, 2),
    )
    assert not molecule.bond_count
    assert _lines(log, 'could not be stored') == [
        'bond: 1 stated bond(s) name an atom that could not be stored, so the bond is not built '
        '(first atoms 0-1)']


# the cross-reader pair


def test_the_same_structure_through_both_readers_is_the_same_molecule():
    """One GLY-ALA dipeptide and one water, read as legacy PDB and as mmCIF, build one molecule.

    The two files reach it by different routes: the legacy file states no bond so every bond is the
    table's, while the mmCIF states its intra-residue bonds in ``_chem_comp_bond`` and they must
    deduplicate against the table.  Neither states the peptide bond and both build it from the
    numbering.  [mutant: `pairs[key] = bond.order` unconditionally in `_stated_bonds`, which doubles
    the mmCIF side's orders and not the legacy side's]
    """
    legacy_log, cif_log = [], []
    legacy = build_molecule(pdb((_DATA / 'pdb_dipeptide.pdb').read_text(encoding='utf-8'))[0], log=legacy_log)
    cif = build_molecule(mmcif((_DATA / 'mmcif_dipeptide.cif').read_text(encoding='utf-8'))[0], log=cif_log)

    assert legacy.atom_count == 11
    assert legacy.bond_count == 9
    assert legacy == cif
    assert legacy_log == cif_log
    # asserting the content and not just the equality is what stops both sides being equal and wrong
    assert legacy_log == []
    # the peptide bond is there in both: two residues and one water, so two components not three
    assert legacy.connected_components_count == 2


# nucleotide end-to-end


def test_a_dinucleotide_builds_the_phosphodiester_link_and_suppresses_the_terminal_op3():
    """A DA-DC dinucleotide comes back bonded, and the mid-chain OP3 suppression fires correctly.

    The nucleotide rows use a distinct ``link_in``/``link_out`` pair (``P`` / ``O3'``).  DC's ``OP3``
    has ``P`` as its sole template neighbour and ``P`` is a built link atom, so it is suppressed; DA is
    5'-terminal with no phosphate at all and the link is built at ``O3'``, so all four absent phosphate
    atoms are reported.  The fixture also writes the aliases ``O3*``, ``O1P`` and ``O2P``, exercising
    ``normalize_atom_name``.  [mutant: dropping the ``one.links_built.add`` call in ``_chain_links``]
    """
    legacy_log = []
    record = pdb((_DATA / 'pdb_dinucleotide.pdb').read_text(encoding='utf-8'))[0]
    mol = build_molecule(record, log=legacy_log)

    assert mol.atom_count == 37       # 18 DA atoms + 19 DC atoms
    assert mol.bond_count == 41       # 20 intra-DA + 20 intra-DC + 1 phosphodiester
    assert mol.connected_components_count == 1

    # the phosphodiester link is O3' of DA (stable id 6) to P of DC (stable id 19), order 1
    bonds = _orders(mol)
    assert bonds.get((6, 19)) == 1, bonds

    # DA (5'-terminal, no phosphate): P, OP1, OP2, OP3 absent and reported
    da_shortfall = _lines(legacy_log, 'DA lack the template atom')
    assert len(da_shortfall) == 1, legacy_log
    assert 'OP1, OP2, OP3, P' in da_shortfall[0]

    # DC (mid-chain): OP3 absent but suppressed, because the link was built at P
    assert not _lines(legacy_log, 'DC lack the template atom'), legacy_log


# the molecule's own log


def test_the_builders_findings_reach_the_molecule_with_nothing_passed_in():
    """`mol.log` is the destination, so a caller who passed no `log=` still has the shortfall report.

    [mutant: drop the `mol.log.absorb('read', own)` in `build_molecule`]
    """
    record = pdb('\n'.join(_residue('ALA', 'A', 1, _ALA[:4])))[0]
    molecule = build_molecule(record)                  # no log= anywhere
    assert _lines(molecule.log, 'ALA lack the template atom'), molecule.log
    assert all(x.stage == 'read' for x in molecule.log), molecule.log


def test_the_records_parse_log_reaches_the_molecule_it_became():
    """A reader finding is about these atoms, so `build_molecule` folds the record's log in too.

    A caller holding the molecule can see the element column was not a symbol without also holding the
    record it came from.  [mutant: absorb `own` only, not `record.log`]
    """
    lines = _residue('ALA', 'A', 1, _ALA)
    lines[0] = _atom(1, 'N', 'ALA', 'A', 1, element=' X')     # column 77-78 is not a symbol
    record = pdb('\n'.join(lines))[0]
    assert _lines(record.log, 'is not an element symbol'), record.log
    molecule = build_molecule(record)
    assert _lines(molecule.log, 'is not an element symbol'), molecule.log


def test_each_model_of_a_multi_model_deck_holds_only_its_own_findings():
    """Model 2's lines are on model 2's molecule and on no other: never pooled.

    Both molecules carry the file-level lines, which the reader gives every record on purpose; the
    per-record line -- which model this is -- is the one that must not cross over.
    [mutant: absorb the caller's flat list instead of this record's]
    """
    deck = ['MODEL        1',
            *_residue('ALA', 'A', 1, _ALA),
            'ENDMDL',
            'MODEL        2',
            *_residue('GLY', 'A', 1, (('N', ' N'), ('CA', ' C')), first=7),
            'ENDMDL']
    records = pdb('\n'.join(deck))
    assert len(records) == 2
    molecules = [build_molecule(record) for record in records]

    assert _lines(molecules[0].log, 'this record holds model 1'), molecules[0].log
    assert not _lines(molecules[0].log, 'this record holds model 2'), molecules[0].log
    assert _lines(molecules[1].log, 'this record holds model 2'), molecules[1].log
    assert not _lines(molecules[1].log, 'this record holds model 1'), molecules[1].log

    # The complete alanine has no shortfall; the two-atom glycine has one, and only there.
    assert not _lines(molecules[0].log, 'lack the template atom'), molecules[0].log
    assert _lines(molecules[1].log, 'GLY lack the template atom'), molecules[1].log


def test_the_caller_list_gets_this_passs_lines_and_the_readers_are_not_repeated_into_it():
    """`log=` receives what this pass found.  The record's own lines are already the reader's caller's,
    so they go to the molecule and are not restated here.

    [mutant: `log.extend(record.log)` beside the absorb]
    """
    lines = _residue('ALA', 'A', 1, _ALA[:4])
    lines[0] = _atom(1, 'N', 'ALA', 'A', 1, element=' X')     # column 77-78 is not a symbol
    record = pdb('\n'.join(lines))[0]
    log: list = []
    molecule = build_molecule(record, log=log)
    assert _lines(log, 'ALA lack the template atom'), log
    assert not _lines(log, 'is not an element symbol'), log
    assert _lines(molecule.log, 'is not an element symbol'), molecule.log


def test_the_mmcif_readers_lines_reach_the_molecule_too():
    """The builder is one pass for both dialects, so an mmCIF record's log lands the same way.

    [mutant: absorb `record.log` in the legacy branch only]
    """
    record = mmcif((_DATA / 'mmcif_ligand_water.cif').read_text(encoding='utf-8'))[0]
    assert record.log                                  # the reader found something to say
    molecule = build_molecule(record)                  # no log= anywhere
    for line in record.log:
        assert _lines(molecule.log, str(line)), (line, molecule.log)


# the ensemble

#: The five heavy atoms of a free glycine, with the element each carries in columns 77-78.
_GLY = (('N', ' N'), ('CA', ' C'), ('C', ' C'), ('O', ' O'), ('OXT', ' O'))


def _model(number, z, atoms=_GLY, first=1):
    """One ``MODEL``/``ENDMDL`` pair holding a free glycine, the whole residue lifted to height ``z``."""
    return [f'{"MODEL":<6}    {number:>4}',
            *(_atom(first + n, name, 'GLY', 'A', 1, x=1.5 * n, y=0.5 * n, z=z, element=element)
              for n, (name, element) in enumerate(atoms)),
            'ENDMDL']


def test_a_list_of_records_collapses_to_one_molecule_with_a_conformer_each():
    """Two models, one molecule, one conformer per model, each carrying its ``MODEL`` number.

    [mutant: build from `records[0]` and ignore the rest]
    """
    records = pdb('\n'.join([*_model(1, 0.25), *_model(2, 0.5, first=6), 'END']))
    assert len(records) == 2
    molecule = build_molecule(records)

    assert len(molecule.conformers) == 2
    assert [c.ext_index for c in molecule.conformers] == [1, 2]
    first, second = molecule.conformers
    assert [round(z, 3) for _, _, z in second.coordinates] == [0.5] * len(molecule)
    assert first.coordinates != second.coordinates
    # Model 0 is what `xyz_of` answers, so the first record is the molecule's own geometry.
    number = molecule.atom_numbers[0]
    assert molecule.xyz_of(number) == first.xyz_of(number)


def test_one_record_still_builds_one_conformer():
    """Every caller today passes one record, so nothing collapses by accident."""
    record = pdb('\n'.join([*_model(1, 0.25), 'END']))[0]
    molecule = build_molecule(record)
    assert len(molecule.conformers) == 1
    assert molecule.conformers[0].ext_index == 1


def test_a_model_with_a_different_atom_set_is_logged_and_skipped():
    """All-or-nothing per model: the short model stores nothing and the first one still lands.

    [mutant: store the atoms that do match and leave the rest at the origin]
    """
    deck = [*_model(1, 0.25), *_model(2, 0.5, atoms=_GLY[:4], first=6), 'END']
    records = pdb('\n'.join(deck))
    log: list = []
    molecule = build_molecule(records, log=log)
    assert len(molecule.conformers) == 1
    assert _lines(log, 'a different atom set'), log


def test_a_model_with_no_number_stores_the_sentinel():
    """A deck with no ``MODEL`` card states no model number, and `ext_index` is None rather than zero."""
    molecule = build_molecule(pdb('\n'.join(_residue('GLY', 'A', 1, _GLY, z=1.))))
    assert molecule.conformers[0].ext_index is None


def test_an_empty_sequence_is_refused():
    """An empty sequence states no molecule to build, which is the caller's error and not a finding."""
    with raises(ValueError, match='at least one record'):
        build_molecule([])
