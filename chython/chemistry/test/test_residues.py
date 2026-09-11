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
"""`tables/residues.tsv` is hand-written connectivity, so it gets an arbiter.

Nothing downstream recomputes it, so a wrong bond order is invisible forever.  Every polymer row is
built from its own atoms and bonds and compared -- container against container, SMILES-string identity
being unsound here -- against a reference SMILES written below independently of the table."""
from pathlib import Path
from re import sub
from subprocess import run
from sys import executable

import pytest

from .._residues import (RESIDUE_KINDS, normalize_atom_name, residue_template, residue_templates)
from .._tables import read_table
from ...core import MoleculeContainer, read_smiles
from ...core._core import element_symbols


#: One reference per polymer row: the free component the Chemical Component Dictionary describes.
#: An amino acid carries both `O` and `OXT`, so it is the amino acid and needs no capping step; a
#: nucleotide carries `OP3`, so it is the 5'-monophosphate.
REFERENCES = {
    'ALA': 'CC(N)C(=O)O',
    'ARG': 'OC(=O)C(N)CCCNC(N)=N',
    'ASN': 'NC(=O)CC(N)C(=O)O',
    'ASP': 'OC(=O)CC(N)C(=O)O',
    'CYS': 'SCC(N)C(=O)O',
    'GLN': 'NC(=O)CCC(N)C(=O)O',
    'GLU': 'OC(=O)CCC(N)C(=O)O',
    'GLY': 'NCC(=O)O',
    'HIS': 'OC(=O)C(N)CC1=CNC=N1',
    'ILE': 'CCC(C)C(N)C(=O)O',
    'LEU': 'CC(C)CC(N)C(=O)O',
    'LYS': 'NCCCCC(N)C(=O)O',
    'MET': 'CSCCC(N)C(=O)O',
    'MSE': 'C[Se]CCC(N)C(=O)O',
    'PHE': 'OC(=O)C(N)CC1=CC=CC=C1',
    'PRO': 'OC(=O)C1CCCN1',
    'SER': 'OCC(N)C(=O)O',
    'THR': 'CC(O)C(N)C(=O)O',
    'TRP': 'OC(=O)C(N)CC1=CNC2=C1C=CC=C2',
    'TYR': 'OC(=O)C(N)CC1=CC=C(O)C=C1',
    'VAL': 'CC(C)C(N)C(=O)O',
    'DA': 'OC1CC(N2C3=C(C(N)=NC=N3)N=C2)OC1COP(=O)(O)O',
    'DC': 'OP(=O)(O)OCC1OC(N2C(=O)N=C(N)C=C2)CC1O',
    'DG': 'OP(=O)(O)OCC1OC(N2C=NC3C(=O)NC(N)=NC=32)CC1O',
    'DT': 'OP(=O)(O)OCC1OC(N2C(=O)NC(=O)C(C)=C2)CC1O',
    'DU': 'OP(=O)(O)OCC1OC(N2C(=O)NC(=O)C=C2)CC1O',
    'A': 'N=1C2=C(N=CN=C2N)N(C=1)C1C(C(C(O1)COP(=O)(O)O)O)O',
    'C': 'OP(=O)(O)OCC1OC(N2C(=O)N=C(N)C=C2)C(O)C1O',
    'G': 'OP(=O)(O)OCC1OC(N2C=NC3C(=O)NC(N)=NC=32)C(O)C1O',
    'U': 'OP(=O)(O)OCC1OC(N2C(=O)NC(=O)C=C2)C(O)C1O',
}

#: The closed scope, spelled out so that widening it is a visible edit rather than a row appearing in
#: a TSV.  Not `SO4`, not `PO4`, not `GOL`, not `HEM`: a polyatomic ligand is `saturate()`'s job.
SCOPE = {
    'amino_acid': {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
                   'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'},
    'nucleotide': {'DA', 'DC', 'DG', 'DT', 'DU', 'A', 'C', 'G', 'U'},
    'water': {'HOH', 'DOD', 'WAT'},
    'ion': {'LI', 'NA', 'K', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'MN', 'MN3', 'FE', 'FE2', 'CO',
            'NI', 'CU', 'CU1', 'ZN', 'CD', 'HG', 'F', 'CL', 'BR', 'IOD'},
}

#: Charge and element for every ion row.  The oxidation state is the only thing distinguishing `FE`
#: from `FE2`, `CU` from `CU1` and `MN` from `MN3`, so without this list a wrong one is invisible.
ION_SPECIES = {
    'LI': ('Li', 1), 'NA': ('Na', 1), 'K': ('K', 1), 'RB': ('Rb', 1), 'CS': ('Cs', 1),
    'MG': ('Mg', 2), 'CA': ('Ca', 2), 'SR': ('Sr', 2), 'BA': ('Ba', 2),
    'MN': ('Mn', 2), 'MN3': ('Mn', 3), 'FE': ('Fe', 3), 'FE2': ('Fe', 2),
    'CO': ('Co', 2), 'NI': ('Ni', 2), 'CU': ('Cu', 2), 'CU1': ('Cu', 1),
    'ZN': ('Zn', 2), 'CD': ('Cd', 2), 'HG': ('Hg', 2),
    'F': ('F', -1), 'CL': ('Cl', -1), 'BR': ('Br', -1), 'IOD': ('I', -1),
}

POLYMERS = sorted(SCOPE['amino_acid'] | SCOPE['nucleotide'])


def build(template) -> MoleculeContainer:
    """The template's own atoms and bonds as a molecule, with implicit hydrogens derived.

    No standardization, no aromatization, no capping -- so a disagreement with a reference is a fact
    about the table.
    """
    molecule = MoleculeContainer()
    with molecule.edit() as edit:
        ids = {name: edit.add_atom(element, charge=charge)
               for name, (element, charge) in template.atoms.items()}
        for a, b, order in template.bonds:
            edit.add_bond(ids[a], ids[b], order)
    molecule.derive_hydrogens()
    return molecule


# --- the arbiter -------------------------------------------------------------------------------- #

@pytest.mark.parametrize('name', POLYMERS)
def test_every_polymer_row_is_the_molecule_its_reference_names(name):
    """Row against hand-written reference, compared as canonical structures."""
    template = residue_template(name)
    assert template is not None, f'{name} is missing from the table'
    built = build(template)
    reference = read_smiles(REFERENCES[name])
    assert built == reference, (
        f'{name} does not match its reference.\n'
        f'  table:     {built.smiles}\n'
        f'  reference: {reference.smiles}\n'
        'The table is what to fix unless the reference can be shown to be the wrong molecule.')


def test_the_arbiter_can_fail():
    """Negative control: a reference deliberately one bond order out must not compare equal.

    Without it the test above is green whenever `==` is doing something other than comparing
    structures.
    """
    built = build(residue_template('PHE'))
    # phenylalanine's ring drawn as cyclohexane: the same atoms, the same connectivity, three orders
    assert built != read_smiles('OC(=O)C(N)CC1CCCCC1')


def test_every_reference_is_flat_and_kekule():
    """No `@` and no lower-case aromatic atom in a reference.

    Both would assert something the table cannot state: it carries no parity, and its rings are Kekule
    because the Chemical Component Dictionary's are.
    """
    for name, smiles in REFERENCES.items():
        assert '@' not in smiles, f'{name}: the table states no parity, so a reference states none'
        # A bracket's contents are exempt: `[Se]` is one element symbol, not an aromatic atom.
        bare = sub(r'\[[^]]*]', '', smiles)
        assert bare == bare.upper(), f'{name}: a lower-case atom is an aromatic one'


# --- the table's own shape ---------------------------------------------------------------------- #

def test_the_scope_is_exactly_what_it_claims():
    templates = residue_templates()
    by_kind = {kind: set() for kind in RESIDUE_KINDS}
    for template in templates.values():
        by_kind[template.kind].add(template.name)
    assert by_kind == SCOPE, (
        'the table\'s scope has drifted.  It is closed on purpose: the Chemical Component '
        'Dictionary has 45000 entries and a table that starts absorbing common ligands has no '
        'natural stopping point.  Widening it is a ruling, not a row.')


def test_there_is_a_reference_for_every_polymer_row_and_nothing_else():
    """The arbiter covers the polymer rows exactly, so a new row cannot arrive uncompared."""
    assert set(REFERENCES) == SCOPE['amino_acid'] | SCOPE['nucleotide']


def test_no_duplicate_component_id():
    """Checked over the file text, because the loader's dict would silently absorb a duplicate."""
    names = [row['name'] for row in read_table('residues.tsv')]
    assert len(names) == len(set(names))


def test_every_bond_and_link_names_a_declared_atom():
    offences = []
    for template in residue_templates().values():
        for a, b, _ in template.bonds:
            for name in (a, b):
                if name not in template.atoms:
                    offences.append(f'{template.name}: bond atom {name}')
        for link in (template.link_in, template.link_out):
            if link is not None and link not in template.atoms:
                offences.append(f'{template.name}: link atom {link}')
    assert not offences, '\n'.join(offences)


def test_no_row_repeats_an_atom_name_or_a_bonded_pair():
    """Over the file text: the loader compiles into dicts and sets, which cannot show a duplicate."""
    offences = []
    for row in read_table('residues.tsv'):
        names = [entry.split(':')[0] for entry in row['atoms'].split(',')]
        if len(names) != len(set(names)):
            offences.append(f'{row["name"]}: repeated atom name')
        pairs = [frozenset(entry.split('-')[:2])
                 for entry in (row['bonds'].split(',') if row['bonds'] else ())]
        if len(pairs) != len(set(pairs)):
            offences.append(f'{row["name"]}: repeated bonded pair')
    assert not offences, '\n'.join(offences)


def test_every_kind_is_one_of_the_four_spellings():
    assert {row['kind'] for row in read_table('residues.tsv')} <= set(RESIDUE_KINDS)


def test_every_element_is_an_element():
    symbols = set(element_symbols()[1:])
    offences = [f'{t.name}:{name}={element}' for t in residue_templates().values()
                for name, (element, _) in t.atoms.items() if element not in symbols]
    assert not offences, '\n'.join(offences)


def test_a_polymer_row_names_both_links_and_a_water_or_ion_names_neither():
    for template in residue_templates().values():
        if template.kind in ('amino_acid', 'nucleotide'):
            assert template.link_in and template.link_out, f'{template.name} links into no chain'
        else:
            assert template.link_in is None and template.link_out is None, template.name


def test_the_polymer_links_are_the_backbone_atoms():
    """The peptide bond is `previous.C -- this.N` and the phosphodiester `previous.O3' -- this.P`.

    The later pass reads these two names and nothing else: a row naming `CA` instead of `N` would still
    load, still bond, and build a chain through the wrong atom.
    """
    for name in sorted(SCOPE['amino_acid']):
        template = residue_template(name)
        assert (template.link_in, template.link_out) == ('N', 'C'), name
    for name in sorted(SCOPE['nucleotide']):
        template = residue_template(name)
        assert (template.link_in, template.link_out) == ('P', "O3'"), name


def test_a_water_or_ion_row_has_no_bonds():
    for template in residue_templates().values():
        if template.kind in ('water', 'ion'):
            assert not template.bonds, template.name


def test_no_bond_in_the_table_is_aromatic():
    """Order 4 is absent by ruling, not by accident, so it is asserted rather than assumed."""
    orders = {order for t in residue_templates().values() for *_, order in t.bonds}
    assert orders <= {1, 2, 3}, f'the table carries orders {sorted(orders)}'


def test_no_polymer_row_carries_a_charge():
    """A PDB file states no protonation state, so a charge on a residue would be invented here.

    Every one of these has a legal neutral valence, so nothing forces the table's hand.
    """
    offences = [f'{t.name}:{name}' for t in residue_templates().values()
                if t.kind in ('amino_acid', 'nucleotide')
                for name, (_, charge) in t.atoms.items() if charge]
    assert not offences, '\n'.join(offences)


def test_every_ion_row_is_one_atom_of_the_species_it_names():
    for name, (element, charge) in ION_SPECIES.items():
        template = residue_template(name)
        assert template is not None and template.kind == 'ion', name
        assert len(template.atoms) == 1, name
        (atom_name, species), = template.atoms.items()
        assert species == (element, charge), f'{name}: {species} is not {(element, charge)}'
        # The CCD spells the oxidation state in the component id and never in the atom name, so an
        # atom of `FE2` is called `FE` and an atom of `IOD` is called `I`.
        assert atom_name in (name, element.upper()), f'{name}: atom named {atom_name}'


def test_water_is_one_oxygen():
    for name in sorted(SCOPE['water']):
        assert residue_template(name).atoms == {'O': ('O', 0)}, name


def test_no_hydrogen_is_listed():
    """Heavy atoms only -- the implicit count comes from the valence rules, per file.

    It is also what makes the terminal cases come out by themselves off one row: a backbone N with two
    heavy neighbours derives one H, the same N at the N-terminus derives two.
    """
    offences = [t.name for t in residue_templates().values()
                if any(element == 'H' for element, _ in t.atoms.values())]
    assert not offences, offences


def test_the_terminal_and_mid_chain_cases_come_off_one_row():
    """Drop `OXT`, add the two peptide-bond neighbours, and the derived counts follow.

    A mid-chain alanine's backbone N holds one hydrogen where the free residue's holds two, and neither
    row nor loader knows which end of a chain it is on.  `OXT` is dropped because a file carries it only
    on the real C-terminus, which is why no bond of it applies mid-chain.
    """
    template = residue_template('ALA')

    free = MoleculeContainer()
    with free.edit() as edit:
        ids = {name: edit.add_atom(element, charge=charge)
               for name, (element, charge) in template.atoms.items()}
        for a, b, order in template.bonds:
            edit.add_bond(ids[a], ids[b], order)
    free.derive_hydrogens()
    assert free.atom(ids['N']).implicit_h == 2

    chained = MoleculeContainer()
    with chained.edit() as edit:
        ids = {name: edit.add_atom(element, charge=charge)
               for name, (element, charge) in template.atoms.items() if name != 'OXT'}
        for a, b, order in template.bonds:
            if 'OXT' not in (a, b):
                edit.add_bond(ids[a], ids[b], order)
        acyl = edit.add_atom('C')                              # the preceding residue's carbonyl
        amide = edit.add_atom('N')                             # the following residue's backbone N
        edit.add_bond(ids['N'], acyl, 1)
        edit.add_bond(ids['C'], amide, 1)
    chained.derive_hydrogens()
    assert chained.atom(ids['N']).implicit_h == 1
    assert chained.atom(ids['C']).implicit_h == 0


# --- the loader --------------------------------------------------------------------------------- #

def _row(*cells) -> dict:
    return dict(zip(('name', 'kind', 'atoms', 'bonds', 'link_in', 'link_out'), cells))


def test_normalize_atom_name():
    assert normalize_atom_name("o3*") == "O3'"          # the 1990s spelling of a ribose oxygen
    assert normalize_atom_name('CA') == 'CA'
    assert normalize_atom_name(' CA ') == 'CA'          # a legacy field is fixed-column and padded
    assert normalize_atom_name("C5'") == "C5'"


def test_the_legacy_phosphate_oxygens_normalize_to_the_modern_spelling():
    """`O1P`/`O2P`/`O3P` are a transposition of `OP1`/`OP2`/`OP3`, so no character rule reaches them.

    The alias applies only for nucleotide rows, because `O1P`, `O2P` and `O3P` are live atom names in
    phosphorylated residues (SEP, TPO, AMP) the table may one day carry.
    """
    assert normalize_atom_name('O1P', kind='nucleotide') == 'OP1'
    assert normalize_atom_name('o2p', kind='nucleotide') == 'OP2'
    assert normalize_atom_name(' O3P ', kind='nucleotide') == 'OP3'
    assert normalize_atom_name('OP1', kind='nucleotide') == 'OP1'   # modern spelling unchanged
    # Without kind the alias is not applied, so the name returns as-is.
    assert normalize_atom_name('O1P') == 'O1P'


def test_a_legacy_spelled_nucleotide_gets_its_whole_phosphate():
    """The defect the alias map exists for, not merely the mapping function.

    A file with the 1990s spellings normalizes `P` and `O3'` fine, so without the map the backbone comes
    out correctly joined with the phosphate hanging off it unbonded.  That looks like it worked, so it is
    pinned as the molecule and not as a string comparison on a name.
    """
    legacy = {'O3P': 'OP3', 'P': 'P', 'O1P': 'OP1', 'O2P': 'OP2', 'O5*': "O5'", 'C5*': "C5'",
              'C4*': "C4'", 'O4*': "O4'", 'C3*': "C3'", 'O3*': "O3'", 'C2*': "C2'", 'C1*': "C1'",
              'N9': 'N9', 'C8': 'C8', 'N7': 'N7', 'C5': 'C5', 'C6': 'C6', 'N6': 'N6', 'N1': 'N1',
              'C2': 'C2', 'N3': 'N3', 'C4': 'C4'}
    template = residue_template('DA')
    # every atom of the row is present in the legacy spelling, and no other -- otherwise a missing
    # atom, not the aliasing, would be what the comparison below detected
    assert set(legacy.values()) == set(template.atoms)

    molecule = MoleculeContainer()
    with molecule.edit() as edit:
        # exactly what the consuming pass does: normalize the file's name, then look the row up
        ids = {}
        for stated in legacy:
            name = normalize_atom_name(stated, kind='nucleotide')
            element, charge = template.atoms[name]
            ids[name] = edit.add_atom(element, charge=charge)
        for a, b, order in template.bonds:
            edit.add_bond(ids[a], ids[b], order)
    molecule.derive_hydrogens()
    assert molecule == read_smiles(REFERENCES['DA'])


def test_no_alias_rewrites_a_name_the_table_already_uses():
    """Every alias entry must be unambiguous in the nucleotide scope it is applied in: a source that is
    a real atom name in a nucleotide row would rewrite that atom away, a target no row uses goes
    nowhere.  Three names are declined by policy, not by these checks: `OW` is GROMACS' `SOL`
    vocabulary, and `OT1`/`OT2` are the CHARMM/XPLOR spellings of `O` and `OXT`, which already live in
    every amino acid row -- aliasing either table-wide would corrupt such a residue.
    """
    from .._residues import _ALIASES

    nucleotide_atoms = {}
    for template in residue_templates().values():
        if template.kind == 'nucleotide':
            for name, (element, _) in template.atoms.items():
                nucleotide_atoms.setdefault(name, set()).add(element)

    for source, target in _ALIASES.items():
        assert source not in nucleotide_atoms, (
            f'{source} is an atom name in a nucleotide row, so aliasing it away loses that atom')
        assert target in nucleotide_atoms, (
            f'{target} is not an atom name any nucleotide row uses; the alias goes nowhere')
        assert len(nucleotide_atoms[target]) == 1, (
            f'{target} names atoms of {sorted(nucleotide_atoms[target])} in different nucleotide rows')

    # Declined by policy: OW (GROMACS/SOL), OT1 and OT2 (CHARMM/XPLOR terminal oxygens)
    assert 'OW' not in _ALIASES
    assert 'OT1' not in _ALIASES
    assert 'OT2' not in _ALIASES


def test_an_alias_does_not_rewrite_a_name_for_its_own_row():
    """Every atom name survives normalisation under its own row's kind unchanged.

    A row whose own atom names include a source of the alias map would be corrupted at load time.  Atom
    names are the table's only join key, so a silent rewrite means a template matching nothing.
    """
    offences = []
    for template in residue_templates().values():
        for atom_name in template.atoms:
            normalised = normalize_atom_name(atom_name, kind=template.kind)
            if normalised != atom_name:
                offences.append(
                    f'{template.name}: {atom_name!r} normalises to {normalised!r} under kind '
                    f'{template.kind!r} -- an atom name in its own row must be stable')
    assert not offences, '\n'.join(offences)


#: Complete atom-name sets for one representative of each `kind`, chosen to cover every naming
#: convention the table uses: ALA for backbone/terminal oxygen/side chain, DA for phosphate/sugar/
#: nucleobase, HOH for water, ZN for a monoatomic ion.  Atom names are the table's only join key, so a
#: silently renamed atom means a bond that never gets built, permanently and invisibly.
PINNED_ATOM_NAMES = {
    'ALA': frozenset({'N', 'CA', 'CB', 'C', 'O', 'OXT'}),
    'DA':  frozenset({'OP3', 'P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'",
                      "C3'", "O3'", "C2'", "C1'",
                      'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'}),
    'HOH': frozenset({'O'}),
    'ZN':  frozenset({'ZN'}),
}


@pytest.mark.parametrize('residue_name', sorted(PINNED_ATOM_NAMES))
def test_the_atom_names_are_exactly_as_pinned(residue_name):
    """Atom names for a representative of each kind, pinned as a frozenset.

    The comparison catches both a rename (`CB` -> `CB1`) and a deletion; either makes the template's
    bonds unreachable for a file using the standard name.
    """
    template = residue_template(residue_name)
    assert template is not None
    assert frozenset(template.atoms) == PINNED_ATOM_NAMES[residue_name], (
        f'{residue_name}: atom-name set differs from the pinned set.  Atom names are the join key '
        'between a file\'s atoms and the template\'s bonds; renaming one silently breaks every '
        'PDB file that carries it.')


def test_an_unknown_component_id_is_a_miss_and_not_a_raise():
    """A ligand or modified base is the common case, so reporting it belongs to the calling pass.

    Raising here would make one unrecognised residue abort a file that read perfectly well.
    """
    assert residue_template('HEM') is None
    assert residue_template('SO4') is None
    assert residue_template('') is None


def test_the_id_is_accepted_in_any_case_and_stripped():
    assert residue_template(' ala ') is residue_template('ALA')
    assert residue_template('hoh') is residue_template('HOH')


def test_the_table_loads_lazily():
    """Importing the package must not read the table, so the check runs in a fresh interpreter.

    In this one the cache is long since populated by the tests above.
    """
    script = ('import chython.chemistry as c\n'
              'from chython.chemistry import _residues\n'
              'assert not _residues._RESIDUES_CACHE, _residues._RESIDUES_CACHE\n'
              'assert _residues.residue_template("ALA") is not None\n'
              'assert _residues._RESIDUES_CACHE\n')
    result = run([executable, '-c', script], capture_output=True, text=True,
                 cwd=str(Path(__file__).resolve().parents[3]))
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('cells,fragment', [
    (('ALA', 'peptide', 'N:N', '', '', ''), 'kind'),
    (('ALA', 'ion', 'NA:Na:1,NA:Na:1', '', '', ''), 'twice'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-CB-1', 'N', 'CA'), 'does not declare'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-CA-1,CA-N-2', 'N', 'CA'), 'bonded twice'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-CA-4', 'N', 'CA'), 'Kekule'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-CA-1', '', 'CA'), 'must name link_in'),
    (('ALA', 'amino_acid', 'N:N:1,CA:C', 'N-CA-1', 'N', 'CA'), 'neutral free component'),
    (('ZN', 'ion', 'ZN:Zn:2', '', 'ZN', ''), 'no neighbour to link to'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-CA', 'N', 'CA'), 'NAME_A-NAME_B-ORDER'),
    (('ALA', 'amino_acid', 'N', '', 'N', 'N'), 'NAME:ELEMENT'),
    (('ALA', 'amino_acid', 'N:N,CA:C', 'N-N-1', 'N', 'CA'), 'to itself'),
])
def test_the_loader_refuses_a_malformed_row(monkeypatch, cells, fragment):
    """Every guard in the loader, exercised: a guard that cannot be shown to fire is a comment."""
    from .. import _residues

    monkeypatch.setattr(_residues, 'read_table', lambda name: [_row(*cells)])
    with pytest.raises(ValueError, match=fragment):
        _residues._compile()


def test_the_loader_accepts_the_rows_it_is_given(monkeypatch):
    """The negative control for the parametrization above: a well-formed pair of rows loads.

    Without it each of those could be raising for a reason unrelated to the guard it names -- a stub
    `read_table` returning the wrong shape, say.
    """
    from .. import _residues

    rows = [_row('ALA', 'amino_acid', 'N:N,CA:C', 'N-CA-1', 'N', 'CA'),
            _row('IOD', 'ion', 'I:I:-1', '', '', '')]
    monkeypatch.setattr(_residues, 'read_table', lambda name: rows)
    out = _residues._compile()
    assert set(out) == {'ALA', 'IOD'}
    assert out['IOD'].atoms == {'I': ('I', -1)}
    assert out['ALA'].bonds == (('N', 'CA', 1),)
    assert (out['ALA'].link_in, out['ALA'].link_out) == ('N', 'CA')


def test_a_duplicate_component_id_is_refused(monkeypatch):
    """An id is the only handle a caller has on a row, so two rows sharing one hides one of them."""
    from .. import _residues

    row = _row('ZN', 'ion', 'ZN:Zn:2', '', '', '')
    monkeypatch.setattr(_residues, 'read_table', lambda name: [row, dict(row)])
    with pytest.raises(ValueError, match='appears twice'):
        _residues._compile()
