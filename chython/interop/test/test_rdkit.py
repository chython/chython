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
"""
The RDKit converter, both directions.  The oracle is `rd0 -> chython 3 -> rd1` judged by InChI twice:
the pair losing nothing, plus `molecule_to_inchi(v3)` read off the arena's own parity bytes, which is
the independent statement a convention error present in both directions cannot cancel.
"""
from csv import DictReader
from pathlib import Path

from pytest import fixture, mark, raises

from chython.core import (H_UNKNOWN, MoleculeContainer as V3Molecule, QueryContainer as V3Query,
                          STEREO_ABS, STEREO_AND, STEREO_OR, WEDGE_UP, molecule_to_inchi)
from chython.exceptions import UnconvertibleType
from .conftest import requires_rdkit
from .._rdkit import from_rdkit, to_rdkit


pytestmark = requires_rdkit

# Corpora: all three optional and all three public.  Two ship inside RDKit; the third is the repo's.
REPO = Path(__file__).resolve().parents[3] / 'test'

# The four stereo unit kinds `_rdkit.py` restates from `chython/core/_stereo.pxi`, each with a
# molecule that can only produce that one, so a renumbering in the core fails here.
KINDS = [('C[C@H](N)O', 0, 'tetrahedral'), ('C/C=C/C', 1, 'cis/trans'),
         ('CC(F)=C=C(F)C', 2, 'allene'), ('Cc1ccccc1-c1ccccc1C', 3, 'atropisomer')]

# Public compounds, one per constitutional feature this converter has to carry.
CONSTITUTION = ['CCO', 'CC(=O)O', 'CC#N', 'C1CC1', 'C1=CC=CC=C1', 'c1ccccc1', 'c1ccncc1',
                'c1cc[nH]c1', 'CN1C=CN=C1', '[13CH4]', '[2H]O[2H]', '[CH3]', '[Na+].[Cl-]',
                'C[N+](C)(C)C.[Cl-]', 'O=S(=O)(O)O', '[O-][N+](=O)c1ccccc1',
                'FC(F)(F)S(=O)(=O)O', 'OC(=O)c1ccccc1O']

TETRAHEDRAL = ['N[C@@H](C)C(=O)O', '[C@@H](N)(C)C(=O)O', 'F[C@](Cl)(Br)I', 'F[C@@](Cl)(Br)I',
               'C[C@H](O)[C@@H](N)CC', 'O[C@H]1CC[C@@H](N)CC1', 'C[S@](=O)c1ccccc1']
CIS_TRANS = ['C/C=C/C', 'C/C=C\\C', 'F/C=C/F', 'F/C=C\\F', 'CC(/C=C/Cl)=C(C)C',
             'O=C(O)/C=C\\C(=O)O', 'C/C(F)=N/O']


def chem():
    """RDKit's `Chem`, quiet.  Imported here so collection without RDKit is a skip, not an error."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog('rdApp.*')
    return Chem


def strip(rd):
    """A copy with no conformers.

    `MolToInchi` derives double-bond geometry from 2D coordinates while `molecule_to_inchi` passes
    none, so a layout-only geometry would give one caller a `/b` layer and not the other.
    """
    Chem = chem()
    out = Chem.Mol(rd)
    out.RemoveAllConformers()
    return out


def inchi_of(rd):
    return chem().MolToInchi(strip(rd))


def layers(text, keys):
    """The named InChI layers of a string, as a dict, so a comparison can name what it compares."""
    return {p[0]: p for p in text.split('/')[1:] if p[0] in keys}


def roundtrip(rd0):
    """`(v3, rd1)` plus the two assertions this file's oracle is made of.

    Only `/t`, `/m` and `/s` are compared against the core: `/h` differs on 55 of the 4991 NCI
    records -- tautomeric N-heterocycles where the core pins what RDKit calls a mobile group.
    """
    v3 = from_rdkit(rd0)
    rd1 = to_rdkit(v3)
    i0, i1 = inchi_of(rd0), molecule_to_inchi(v3)
    assert i0 == inchi_of(rd1), 'the pair lost something'
    assert layers(i0, 'tms') == layers(i1, 'tms'), 'the core disagrees about the configuration'
    return v3, rd1


def state(mol):
    """Every atom property this converter claims to carry, in the molecule's own order.

    A second oracle: InChI has no radical layer and normalizes charges.
    """
    return [(a.element, a.charge, a.isotope, a.is_radical, a.implicit_h, a.map_number)
            for a in mol.atoms()]


def test_atom_order_is_the_molecules_own_v3():
    """RDKit index order is `atoms()` order and nothing else.

    `depict/layout/molecule.py:53` reads coordinates back by position.  The molecule is built so its
    iteration order is neither sorted nor its stable-id order: ids 1, 3, 4, 5.
    """
    mol = V3Molecule()
    with mol.edit():
        ids = [mol.add_atom(6) for _ in range(4)]
        for n, m in zip(ids, ids[1:]):
            mol.add_bond(n, m, 1)
    with mol.edit():
        mol.delete_atom(ids[1])
    with mol.edit():
        extra = mol.add_atom(8)
        mol.add_bond(ids[2], extra, 1)

    order = [a.n for a in mol.atoms()]
    assert order == [1, 3, 4, 5]
    # `keep_numbers` is how the ids become visible on the RDKit side at all; the position they land in
    # is the claim under test, and it is the same whatever is written in the map field.
    assert [a.GetAtomMapNum() for a in to_rdkit(mol, keep_numbers=True).GetAtoms()] == order


@mark.parametrize('text', ['c1ccccc1', 'c1ccncc1', 'c1cc[nH]c1', 'c1ccc2ccccc2c1'])
def test_aromatic_bonds_are_exported_as_aromatic(text):
    """A stored order-4 bond goes out aromatic, on the bond and on both its atoms."""
    Chem = chem()
    mol = from_rdkit(Chem.MolFromSmiles(text))
    aromatic = {b.n for b in mol.bonds() if b.order == 4} | {b.m for b in mol.bonds() if b.order == 4}
    assert aromatic, 'the fixture is not aromatic, so this test would prove nothing'

    rd = to_rdkit(mol, keep_numbers=True)  # so a map number in the export names a chython atom
    for b in rd.GetBonds():
        n, m = b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()
        if mol.order_of(n, m) == 4:
            assert b.GetBondType() == Chem.BondType.AROMATIC and b.GetIsAromatic()
        else:
            assert b.GetBondType() != Chem.BondType.AROMATIC
    assert {a.GetAtomMapNum() for a in rd.GetAtoms() if a.GetIsAromatic()} == aromatic


def test_kekule_input_stays_kekule():
    """Nothing here aromatizes what RDKit handed over kekulized."""
    Chem = chem()
    rd = Chem.MolFromSmiles('c1ccccc1')
    Chem.Kekulize(rd, clearAromaticFlags=True)
    assert sorted(b.order for b in from_rdkit(rd).bonds()) == [1, 1, 1, 2, 2, 2]


def test_aromaticity_is_read_from_the_bond_type_not_the_flag():
    """The import direction reads `GetBondType()` and never `GetIsAromatic()`.

    The two are separable: this molecule has aromatic bond types with every aromatic flag cleared, so
    reading the flag would store six single bonds.
    """
    Chem = chem()
    rd = Chem.RWMol(Chem.MolFromSmiles('c1ccccc1'))
    for a in rd.GetAtoms():
        a.SetIsAromatic(False)
    for b in rd.GetBonds():
        b.SetIsAromatic(False)
    assert all(b.GetBondType() == Chem.BondType.AROMATIC for b in rd.GetBonds())
    assert sorted(b.order for b in from_rdkit(rd).bonds()) == [4] * 6


def test_a_loss_is_reported_and_not_raised():
    """A molecule RDKit cannot hold in full still converts, and says what it dropped."""
    mol = V3Molecule()
    with mol.edit():
        centre = mol.add_atom(6, implicit_h=0)
        for _ in range(5):  # a valence RDKit refuses
            mol.add_bond(centre, mol.add_atom(6, implicit_h=3), 1)

    log = []
    rd = to_rdkit(mol, log=log)
    assert rd.GetNumAtoms() == 6, 'the molecule is the caller\'s, garbage or not'
    assert any('sanitization failed' in x for x in log)
    to_rdkit(mol)  # and log is optional: no list, no raise, same molecule


@mark.parametrize('value', ['CCO', 42, None, b'CCO', object()])
def test_the_import_direction_refuses_what_it_cannot_read(value):
    """`UnconvertibleType`, naming the type, not an `AttributeError`."""
    with raises(UnconvertibleType):
        from_rdkit(value)


def test_the_export_direction_refuses_a_query():
    """A query is neither a molecule nor a reaction and there is no RDKit form to give it.

    It passes `interop.is_container` and so arrives at the exporter rather than at the importer; the
    export side is what refuses it.  A REACTION IS NOT IN THIS LIST ANY MORE: it converts as a
    `ChemicalReaction`, which is what `test_a_reaction_round_trips_through_a_chemical_reaction` asserts.
    """
    from chython.core import read_smarts

    for bad in (read_smarts('[C;D2]'), V3Query()):
        with raises(UnconvertibleType):
            to_rdkit(bad)


def test_the_public_path_reaches_both_directions():
    """`interop.rdkit` is the name callers use, so it is asserted rather than assumed."""
    from chython.core import read_smiles
    from chython.interop import rdkit

    Chem = chem()
    assert isinstance(rdkit(read_smiles('CCO')), Chem.Mol)
    assert isinstance(rdkit(Chem.MolFromSmiles('CCO')), V3Molecule)


def test_the_container_methods_answer_what_the_function_does():
    """`mol.to_rdkit()` and `rxn.to_rdkit()` against a real RDKit, keywords and all.

    `chython/core/test/test_interop_injection.py` proves the wiring with stubs; this is the same two
    methods with the toolkit actually present, which is what a caller in a notebook types.
    """
    from chython.core import read_smiles

    Chem = chem()
    mol = read_smiles('CCO')
    assert Chem.MolToSmiles(mol.to_rdkit()) == Chem.MolToSmiles(to_rdkit(mol))
    assert ([a.GetAtomMapNum() for a in mol.to_rdkit(keep_numbers=True).GetAtoms()]
            == [a.n for a in mol.atoms()])

    rxn = read_smiles('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]')
    rd = rxn.to_rdkit()
    assert Chem.rdChemReactions.ReactionToSmiles(rd) == \
           Chem.rdChemReactions.ReactionToSmiles(to_rdkit(rxn))
    assert ':1' in Chem.rdChemReactions.ReactionToSmiles(rd)


@mark.parametrize('text', CONSTITUTION)
def test_constitution_round_trip(text):
    """Elements, isotopes, charges, radicals, bond orders and hydrogen counts, both ways.

    The stored-state comparison is the second oracle: InChI has no radical layer at all.
    """
    Chem = chem()
    rd0 = Chem.MolFromSmiles(text)
    log = []
    v3 = from_rdkit(rd0, log=log)
    assert log == []
    again, _ = roundtrip(rd0)
    assert state(v3) == state(again)
    # canonical SMILES as well as InChI, because InChI normalizes and SMILES does not.  No keyword:
    # an unmapped molecule exports with an empty map field, so the two SMILES are comparable as they
    # come -- which is the default this asserts.
    assert Chem.MolToSmiles(to_rdkit(v3)) == Chem.MolToSmiles(rd0)


def test_radicals_survive_the_round_trip():
    """The one constitutional feature InChI is blind to, asserted on its own."""
    Chem = chem()
    v3, rd1 = roundtrip(Chem.MolFromSmiles('[CH3]'))
    assert [a.is_radical for a in v3.atoms()] == [True]
    assert [a.GetNumRadicalElectrons() for a in rd1.GetAtoms()] == [1]


def test_dative_bonds_point_from_the_ligand_to_the_metal():
    """Order 8 becomes `BondType.DATIVE`, and RDKit's dative bond is directed while chython's is not.

    `_MAIN_GROUP` is what reconstructs the direction; without it the bond goes out Fe -> N.
    """
    Chem = chem()
    mol = V3Molecule()
    with mol.edit():
        n = mol.add_atom(7, implicit_h=3)
        fe = mol.add_atom(26, implicit_h=0)
        mol.add_bond(fe, n, 8)  # deliberately metal-first, so the export has to swap it

    bond = to_rdkit(mol).GetBondWithIdx(0)
    assert bond.GetBondType() == Chem.BondType.DATIVE
    assert (bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()) == ('N', 'Fe')
    assert [b.order for b in from_rdkit(to_rdkit(mol)).bonds()] == [8]


def test_an_unmodelled_bond_type_becomes_order_eight_and_is_named():
    """A bond chython has no order for becomes order 8, not a single bond the source never stated."""
    Chem = chem()
    rd = Chem.RWMol(Chem.MolFromSmiles('CC'))
    rd.GetBondWithIdx(0).SetBondType(Chem.BondType.QUADRUPLE)
    log = []
    assert [b.order for b in from_rdkit(rd, log=log).bonds()] == [8]
    assert any('QUADRUPLE' in x for x in log)


def test_h_unknown_is_reported_and_never_becomes_zero():
    """`H_UNKNOWN` is reported, and never becomes a stated zero, which is a different molecule.

    The atom goes out with no count at all: `SetNoImplicit` is not called, so RDKit perceives one.
    """
    mol = V3Molecule()
    with mol.edit():
        a = mol.add_atom(6, implicit_h=H_UNKNOWN)
        mol.add_bond(a, mol.add_atom(6, implicit_h=3), 1)

    log = []
    rd = to_rdkit(mol, log=log)
    assert any('H_UNKNOWN' in x for x in log)
    assert [(a.GetNumExplicitHs(), a.GetNoImplicit()) for a in rd.GetAtoms()] == \
           [(0, False), (3, True)]


def test_hydrogen_counts_are_stated_by_default_and_perceived_on_request():
    """`keep_hydrogens` decides who owns the count.

    True states it with `SetNoImplicit` so RDKit cannot add to it; re-perception loses a hydrogen.
    """
    Chem = chem()
    mol = from_rdkit(Chem.MolFromSmiles('c1cc[nH]c1'))
    assert all(a.GetNoImplicit() for a in to_rdkit(mol).GetAtoms())
    assert not any(a.GetNoImplicit() for a in to_rdkit(mol, keep_hydrogens=False).GetAtoms())


def test_keep_mapping_is_the_atom_atom_mapping_and_nothing_else():
    """`keep_mapping` carries `map_number`, so an unmapped molecule exports with an empty map field.

    The half a caller reads back: an export of a mapped record and nothing added to an unmapped one.
    """
    Chem = chem()
    mol = from_rdkit(Chem.MolFromSmiles('CCO'))
    assert [a.map_number for a in mol.atoms()] == [0, 0, 0], 'the fixture is unexpectedly mapped'
    assert [a.GetAtomMapNum() for a in to_rdkit(mol).GetAtoms()] == [0, 0, 0]

    with mol.edit():
        mol.set_map_number(next(iter(mol.atoms())).n, 42)
    assert 42 in [a.GetAtomMapNum() for a in to_rdkit(mol).GetAtoms()]
    assert [a.GetAtomMapNum() for a in to_rdkit(mol, keep_mapping=False).GetAtoms()] == [0, 0, 0]


def test_keep_numbers_writes_the_stable_id_and_says_what_it_displaced():
    """`keep_numbers` puts chython's own atom ids in the field instead: a label, not a mapping.

    RDKit has one integer per atom, so a `map_number` that is not the atom number cannot also fit; the
    conflict is logged rather than guessed at, and only when there is one to report.
    """
    Chem = chem()
    mol = from_rdkit(Chem.MolFromSmiles('CCO'))
    assert ([a.GetAtomMapNum() for a in to_rdkit(mol, keep_numbers=True).GetAtoms()]
            == [a.n for a in mol.atoms()])

    log = []
    to_rdkit(mol, keep_numbers=True, log=log)
    assert not any('map number' in x for x in log), 'nothing was displaced, so nothing to report'

    with mol.edit():
        mol.set_map_number(next(iter(mol.atoms())).n, 42)
    log = []
    assert 42 not in [a.GetAtomMapNum() for a in to_rdkit(mol, keep_numbers=True, log=log).GetAtoms()]
    assert any('map number' in x for x in log)
    # and no complaint when the mapping was not asked for: nothing was competing for the field
    log = []
    to_rdkit(mol, keep_numbers=True, keep_mapping=False, log=log)
    assert not any('map number' in x for x in log)


def test_a_reaction_round_trips_through_a_chemical_reaction():
    """The three sides keep their sides and the atom-atom mapping survives both directions."""
    from chython.core import ReactionContainer, read_smiles

    Chem = chem()
    rxn = read_smiles('[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]>'
                      'O.[Na+].[OH-]>'
                      '[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]')
    assert isinstance(rxn, ReactionContainer)

    rd = to_rdkit(rxn)
    assert (rd.GetNumReactantTemplates(), rd.GetNumAgentTemplates(),
            rd.GetNumProductTemplates()) == (2, 3, 1)

    back = from_rdkit(rd)
    assert isinstance(back, ReactionContainer)
    assert [len(m) for m in back.reactants] == [len(m) for m in rxn.reactants]
    assert [len(m) for m in back.agents] == [len(m) for m in rxn.agents]
    assert ([[a.map_number for a in m.atoms()] for m in back.molecules()]
            == [[a.map_number for a in m.atoms()] for m in rxn.molecules()])
    # and RDKit's own writer prints the mapping, which is why the reaction goes through this type
    assert ':1' in Chem.rdChemReactions.ReactionToSmiles(rd)


def test_an_unmapped_reaction_exports_unmapped():
    """`keep_mapping` means the same thing per molecule inside a reaction as it does outside one."""
    from chython.core import read_smiles

    rxn = read_smiles('CC(=O)O.CCO>>CC(=O)OCC')
    rd = to_rdkit(rxn)
    assert not any(a.GetAtomMapNum() for m in rd.GetReactants() for a in m.GetAtoms())
    assert not any(a.GetAtomMapNum() for m in rd.GetProducts() for a in m.GetAtoms())


def test_a_reaction_shares_one_log_with_its_molecules():
    """One list reports the whole record: a molecule's loss inside a reaction is the reaction's loss."""
    from chython.core import ReactionContainer, read_smiles

    unknown = V3Molecule()
    with unknown.edit():
        a = unknown.add_atom(6, implicit_h=H_UNKNOWN)
        unknown.add_bond(a, unknown.add_atom(6, implicit_h=3), 1)

    log = []
    to_rdkit(ReactionContainer([unknown], [read_smiles('CC')]), log=log)
    assert any('H_UNKNOWN' in x for x in log)


def test_out_of_range_values_are_clamped_and_named():
    """Charges beyond -4..8, radical counts above one and map numbers above 9999 do not fit.

    Each is clamped or dropped and each says so; the fourth line is this molecule's absent valence
    cache.
    """
    Chem = chem()
    rd = Chem.RWMol()
    a = Chem.Atom(6)
    a.SetFormalCharge(9)
    a.SetNumRadicalElectrons(2)
    a.SetAtomMapNum(10001)
    rd.AddAtom(a)

    log = []
    mol = from_rdkit(rd, log=log)
    atom = next(iter(mol.atoms()))
    assert (atom.charge, atom.is_radical, atom.map_number) == (8, True, 0)
    assert len(log) == 4 and any('valence cache' in x for x in log)
    assert any('clamped' in x and 'charge' in x for x in log)
    assert any('radical' in x for x in log)
    assert any('map number' in x for x in log)


@mark.parametrize('text,kind,name', KINDS, ids=[x[2] for x in KINDS])
def test_stereo_unit_kinds_are_what_the_converter_thinks(text, kind, name):
    """The four kind numbers `_rdkit.py` restates from the core, one molecule per kind.

    The core does not publish these on its Python surface, so a renumbering would be silent: a
    cis/trans unit read as a tetrahedron writes a chiral tag onto a double bond.
    """
    mol = from_rdkit(chem().MolFromSmiles(text))
    assert {u['kind'] for u in mol.stereo_units() if u['stereogenic']} == {kind}, name


@mark.parametrize('text', TETRAHEDRAL)
def test_tetrahedral_round_trip(text):
    """Both frame calibrations at once, on every stereocentre.

    RDKit lists an atom's neighbours in bond order with the implicit hydrogen appended last, and
    `CHI_TETRAHEDRAL_CCW` is SMILES `@` over that list.
    """
    v3, _ = roundtrip(chem().MolFromSmiles(text))
    assert any(u['parity'] and u['kind'] == 0 for u in v3.stereo_units())


def test_the_hydrogen_is_last_in_rdkits_neighbour_order():
    """The hydrogen-last calibration, isolated.

    These two strings are one molecule with the implicit hydrogen at opposite ends of the SMILES
    ordering, and RDKit reports opposite chiral tags; only a hydrogen-last rule reconciles that.
    """
    Chem = chem()
    assert inchi_of(Chem.MolFromSmiles('N[C@@H](C)C(=O)O')) != \
           inchi_of(Chem.MolFromSmiles('[C@@H](N)(C)C(=O)O')), 'the fixture stopped separating them'
    for text in ('N[C@@H](C)C(=O)O', '[C@@H](N)(C)C(=O)O'):
        rd0 = Chem.MolFromSmiles(text)
        assert molecule_to_inchi(from_rdkit(rd0)) == inchi_of(rd0), text


@mark.parametrize('text', CIS_TRANS)
def test_cis_trans_round_trip(text):
    """`test/stereo.sdf` carries no double-bond geometry, so these are hand-written and public.

    They are the only anchor for the STEREOZ/STEREOE polarity, which `roundtrip` does not cover, so
    `/b` is asserted here.
    """
    rd0 = chem().MolFromSmiles(text)
    v3, _ = roundtrip(rd0)
    i0 = inchi_of(rd0)
    assert layers(i0, 'b'), 'the fixture carries no geometry, so this would prove nothing'
    assert layers(i0, 'b') == layers(molecule_to_inchi(v3), 'b'), 'the core disagrees about /b'


@mark.parametrize('text,name', [('CC(F)=C=C(F)C', 'allene/cumulene'),
                                ('Cc1ccccc1-c1ccccc1C', 'atropisomer')])
def test_kinds_rdkit_cannot_hold_are_dropped_and_said(text, name):
    """RDKit has no form for an allene or an atropisomer, so exporting one is a stated loss.

    The parity is set by hand because no RDKit input can produce one, which is the point.
    """
    mol = from_rdkit(chem().MolFromSmiles(text))
    unit = next(u for u in mol.stereo_units() if u['stereogenic'] and u['kind'] in (2, 3))
    with mol.edit():
        mol.set_parity(unit['anchor'], 1)

    log = []
    to_rdkit(mol, log=log)
    assert [str(x) for x in log] == [f'1 {name} configuration(s) dropped: RDKit has no form for them']


def test_a_wedge_with_no_parity_is_reported():
    """A drawn centre nobody derived a parity for reaches RDKit as no configuration, and is logged.

    RDKit takes tags, not wedges, and the core stores the wedge until the derivation lands.
    """
    mol = from_rdkit(chem().MolFromSmiles('CC(N)O'))
    unit = next(iter(mol.stereogenic_units()))
    with mol.edit():
        mol.set_wedge(unit['anchor'], unit['refs'][0], WEDGE_UP)

    log = []
    to_rdkit(mol, log=log)
    assert any('wedge' in x for x in log)


def test_unknown_geometry_and_non_tetrahedral_tags_are_dropped_and_said():
    """RDKit's "either" double bond and its non-tetrahedral chiral tags are dropped and said.

    A silently ignored `STEREOANY` reads back as "geometry not specified", a weaker statement.
    """
    Chem = chem()
    rd = Chem.RWMol(Chem.MolFromSmiles('CC=CC'))
    rd.GetBondWithIdx(1).SetStereo(Chem.BondStereo.STEREOANY)
    log = []
    from_rdkit(rd, log=log)
    assert any('unknown geometry' in x for x in log)

    rd2 = Chem.RWMol(Chem.MolFromSmiles('F[Pt](F)(F)F'))
    rd2.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
    log2 = []
    from_rdkit(rd2, log=log2)
    assert any('not tetrahedral' in x for x in log2)


def test_stereo_groups_round_trip_and_keep_their_ids():
    """AND and OR groups both ways, with the group ids preserved.

    RDKit renumbers groups from one on write unless `SetWriteId` says otherwise.
    """
    Chem = chem()
    rd0 = Chem.MolFromSmiles('C[C@H](O)[C@@H](N)CC |o2:3,&1:1|')
    log = []
    v3 = from_rdkit(rd0, log=log)
    assert log == [], 'nothing had to be renumbered: both ids fit'
    groups = v3.stereo_groups()
    assert {k: len(m) for k, m in groups.items()} == {(STEREO_AND, 1): 1, (STEREO_OR, 2): 1}

    rd1 = to_rdkit(v3)
    assert from_rdkit(rd1).stereo_groups() == groups
    cx = Chem.MolToCXSmiles(rd1)
    assert '|o2:' in cx and '&1:' in cx


def test_the_absolute_group_survives_the_import():
    """An ABS group is imported, not dropped: "known rather than a mixture" is a statement."""
    Chem = chem()
    from rdkit.Chem import CreateStereoGroup, StereoGroupType

    rd = Chem.RWMol(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)CC'))
    rd.SetStereoGroups([CreateStereoGroup(StereoGroupType.STEREO_ABSOLUTE, rd, [1], [])])
    assert list(from_rdkit(rd).stereo_groups()) == [(STEREO_ABS, 0)]


def test_out_of_range_group_ids_are_renumbered_and_said():
    """chython stores 1..63 per kind and RDKit's ids are free integers, so some are renumbered.

    Two groups sharing one out-of-range id must not collapse into one.
    """
    Chem = chem()
    from rdkit.Chem import CreateStereoGroup, StereoGroupType

    rd = Chem.RWMol(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)CC'))
    rd.SetStereoGroups([CreateStereoGroup(StereoGroupType.STEREO_AND, rd, [1], [], 999),
                        CreateStereoGroup(StereoGroupType.STEREO_AND, rd, [3], [], 999)])
    log = []
    groups = from_rdkit(rd, log=log).stereo_groups()
    assert sorted(groups) == [(STEREO_AND, 1), (STEREO_AND, 2)]
    assert [str(x) for x in log] == ['2 stereo group id(s) renumbered: chython stores 1..63 per kind']


def test_an_atom_in_two_groups_keeps_the_first():
    """RDKit allows one atom in an AND group and an OR group at once; chython stores one mark."""
    Chem = chem()
    from rdkit.Chem import CreateStereoGroup, StereoGroupType

    rd = Chem.RWMol(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)CC'))
    rd.SetStereoGroups([CreateStereoGroup(StereoGroupType.STEREO_AND, rd, [1], [], 1),
                        CreateStereoGroup(StereoGroupType.STEREO_OR, rd, [1], [], 1)])
    assert list(from_rdkit(rd).stereo_groups()) == [(STEREO_AND, 1)]


def test_absolute_keyword_names_the_unclaimed_centres():
    """`absolute=True` adds an ABS group over the stereocentres no AND/OR group claims.

    Only over centres whose configuration was written: RDKit drops a group with no chiral tag.
    """
    Chem = chem()
    from rdkit.Chem import StereoGroupType

    mol = from_rdkit(Chem.MolFromSmiles('C[C@H](O)[C@@H](N)CC |&1:1|'))
    rd = to_rdkit(mol, absolute=True)
    assert {g.GetGroupType(): len(g.GetAtoms()) for g in rd.GetStereoGroups()} == \
           {StereoGroupType.STEREO_AND: 1, StereoGroupType.STEREO_ABSOLUTE: 1}
    assert not to_rdkit(from_rdkit(Chem.MolFromSmiles('CCO')), absolute=True).GetStereoGroups()


def test_coordinates_round_trip_as_a_2d_conformer():
    """The layout round-trips as a 2D conformer; `keep_coordinates=None` exports it when there is one.

    `Set3D(False)` is not cosmetic: RDKit's depiction and its MDL writer both read the flag.
    """
    mol = V3Molecule()
    with mol.edit():
        a = mol.add_atom(6)
        b = mol.add_atom(8)
        mol.add_bond(a, b, 1)
    with mol.edit():
        mol.set_xy(a, 1.5, -2.5)
        mol.set_xy(b, 3., 0.)

    rd = to_rdkit(mol)
    assert rd.GetNumConformers() == 1 and not rd.GetConformer().Is3D()
    assert [(x.x, x.y) for x in from_rdkit(rd).atoms()] == [(1.5, -2.5), (3., 0.)]
    assert to_rdkit(mol, keep_coordinates=False).GetNumConformers() == 0

    flat = V3Molecule()
    with flat.edit():
        flat.add_atom(6)
    assert to_rdkit(flat).GetNumConformers() == 0, 'no layout, no conformer'
    assert to_rdkit(flat, keep_coordinates=True).GetNumConformers() == 1


POSITIONS = [(0., 0., 0.), (1.5, 0., .3), (3., 0., 0.)]


def _with_conformer(smiles, positions, solid):
    rd = chem().MolFromSmiles(smiles)
    conf = chem().Conformer(rd.GetNumAtoms())
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(solid)
    rd.AddConformer(conf, assignId=True)
    return rd


def test_an_imported_3d_conformer_is_stored_and_is_not_a_layout():
    """An imported 3D conformer is stored in `SEG_CONFORMERS` and is not a loss, but is not a layout.

    The xy of a 3D conformer is a projection, and using one as a drawing would invent a depiction the
    source never had -- so `has_3d` is True while `has_coordinates` stays False.
    """
    log = []
    back = from_rdkit(_with_conformer('CCO', POSITIONS, True), log=log)
    assert back.has_3d is True
    assert [back.xyz_of(n) for n in back] == POSITIONS
    assert not back.has_coordinates, 'a projection is not a layout'
    assert not any('conformer' in x for x in log), f'nothing was lost, so nothing is said: {log}'


def test_a_2d_and_a_3d_conformer_fill_one_segment_each():
    """Two facts, two segments, and neither overwrites the other -- the whole point of the split."""
    rd = _with_conformer('CCO', POSITIONS, True)
    flat = chem().Conformer(rd.GetNumAtoms())
    for i, xy in enumerate([(7., 8., 0.), (9., 10., 0.), (11., 12., 0.)]):
        flat.SetAtomPosition(i, xy)
    flat.Set3D(False)
    rd.AddConformer(flat, assignId=True)

    back = from_rdkit(rd)
    assert back.has_3d and back.has_coordinates
    assert [back.xyz_of(n) for n in back] == POSITIONS
    assert [back.xy_of(n) for n in back] == [(7., 8.), (9., 10.), (11., 12.)]


def test_every_3d_conformer_becomes_a_model_and_nothing_is_a_loss():
    """An ensemble's members are models, in the source's order, each keeping its RDKit conformer id.

    [mutant: keep the first 3D conformer and log the rest as dropped]
    """
    rd = _with_conformer('CCO', POSITIONS, True)
    for shift in (1., 2.):
        conf = chem().Conformer(rd.GetNumAtoms())
        for i, (x, y, z) in enumerate(POSITIONS):
            conf.SetAtomPosition(i, (x + shift, y, z))
        conf.Set3D(True)
        rd.AddConformer(conf, assignId=True)

    log = []
    back = from_rdkit(rd, log=log)
    assert len(back.conformers) == 3
    assert [back.xyz_of(n) for n in back] == POSITIONS, 'model 0 is the first conformer'
    assert back.conformer(2).coordinates == [(x + 2., y, z) for x, y, z in POSITIONS]
    assert [c.ext_index for c in back.conformers] == [c.GetId() for c in rd.GetConformers()]
    assert not any('conformer' in x for x in log), f'nothing was lost, so nothing is said: {log}'


def test_every_model_goes_back_out_as_a_conformer_after_the_layout():
    """Both directions: N models in, N conformers out, and the layout still holds id 0.

    [mutant: export `mol.conformer(0)` alone]
    """
    mol = from_rdkit(_with_conformer('CCO', POSITIONS, False))     # a layout and no geometry
    numbers = mol.atom_numbers
    with mol.edit():
        for k in range(3):
            model = mol.add_conformer(ext_index=k)
            for i, n in enumerate(numbers):
                mol.set_xyz(n, float(k), float(i), 0., model=model)

    out = to_rdkit(mol)
    assert out.GetNumConformers() == 4, 'the layout plus three models'
    assert out.GetConformer(0).Is3D() is False, 'the layout keeps id 0'
    assert [tuple(out.GetConformer(3).GetAtomPosition(i)) for i in range(3)] == [(2., 0., 0.),
                                                                                (2., 1., 0.),
                                                                                (2., 2., 0.)]
    back = from_rdkit(out)
    assert len(back.conformers) == 3
    assert back.conformer(2).xyz_of(back.atom_numbers[1]) == (2., 1., 0.)


def test_the_geometry_is_exported_as_a_second_conformer_and_the_layout_stays_first():
    """`GetConformer()` returns id 0, so the layout must keep that slot.

    A depiction engine calling `GetConformer()` on a 3D molecule must still get the drawing.
    """
    rd = _with_conformer('CCO', POSITIONS, True)
    flat = chem().Conformer(rd.GetNumAtoms())
    for i, xy in enumerate([(7., 8., 0.), (9., 10., 0.), (11., 12., 0.)]):
        flat.SetAtomPosition(i, xy)
    flat.Set3D(False)
    rd.AddConformer(flat, assignId=True)
    mol = from_rdkit(rd)

    out = to_rdkit(mol)
    assert out.GetNumConformers() == 2
    assert out.GetConformer(0).Is3D() is False
    assert out.GetConformer(1).Is3D() is True
    assert [tuple(out.GetConformer(1).GetAtomPosition(i)) for i in range(3)] == POSITIONS


def test_a_geometry_only_molecule_exports_its_geometry_even_with_the_layout_suppressed():
    """`keep_coordinates` is about a depiction and must not gate the geometry.

    A caller suppressing the layout has said nothing about the structure's coordinates.
    """
    mol = from_rdkit(_with_conformer('CCO', POSITIONS, True))
    assert mol.has_3d and not mol.has_coordinates
    out = to_rdkit(mol, keep_coordinates=False)
    assert out.GetNumConformers() == 1
    assert out.GetConformer().Is3D() is True
    assert [tuple(out.GetConformer().GetAtomPosition(i)) for i in range(3)] == POSITIONS


def _rdkit_data():
    """RDKit's own data directory, which is where two of the three corpora live."""
    import rdkit

    return Path(rdkit.__file__).resolve().parent


@fixture(scope='module')
def stereo_sdf():
    """`test/stereo.sdf` as RDKit molecules: 300 records dense in tetrahedral stereo."""
    path = REPO / 'stereo.sdf'
    if not path.is_file():
        from pytest import skip

        skip(f'the repo stereo corpus is not present (looked for {path})')
    from chython.formats.ctfile import SDFRead

    with SDFRead(str(path)) as f:
        return [to_rdkit(m) for m in f]


@fixture(scope='module')
def nci():
    """4999 public NCI records shipped inside RDKit.  No stereo, so this is a constitution sweep."""
    path = _rdkit_data() / 'Data' / 'NCI' / 'first_5K.smi'
    if not path.is_file():
        from pytest import skip

        skip(f'the RDKit NCI corpus is not present (looked for {path})')
    Chem = chem()
    out = []
    # `encoding='utf-8'` here and on the filter table below: RDKit ships both as UTF-8 (the table's
    # notes hold a `‐`), and `open` without it asks the locale -- cp1252 on the Windows runner, which
    # cannot decode a continuation byte and turns a corpus sweep into an error.
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.split():
                rd = Chem.MolFromSmiles(line.split()[0])
                if rd is not None:  # 8 records RDKit itself will not read; nothing to compare
                    out.append(rd)
    assert len(out) > 4900, f'the corpus shrank: {len(out)} records'
    return out


@fixture(scope='module')
def filter_examples():
    """A substructure-filter table in RDKit's contrib tree, read for its PubChem examples.

    1826 unique public structures, 263 with configured double-bond geometry -- the one thing the
    repo's own `test/` has none of.
    """
    path = (_rdkit_data() / 'Contrib' / 'NIBRSubstructureFilters' /
            'SubstructureFilter_HitTriaging_wPubChemExamples.csv')
    if not path.is_file():
        from pytest import skip

        skip(f'the RDKit substructure-filter table is not present (looked for {path})')
    Chem = chem()
    texts = set()
    with path.open(encoding='utf-8') as f:
        for row in DictReader(f):
            for column in ('EX1', 'EX2', 'EX3', 'EX4', 'EX5'):
                text = (row.get(column) or '').strip()
                if text:
                    texts.add(text)
    out = [Chem.MolFromSmiles(x) for x in sorted(texts)]
    out = [x for x in out if x is not None]
    assert len(out) > 1800, f'the table shrank: {len(out)} records'
    return out


@mark.parametrize('corpus', ['stereo_sdf', 'nci', 'filter_examples'])
def test_corpus_round_trips(corpus, request):
    """`rd0 -> chython 3 -> rd1` loses nothing, over 7117 records from three public corpora.

    The pair assertion, so it is blind to a convention error present in both directions;
    `test_corpus_configurations_match_the_core` is the half that is not.
    """
    bad = []
    for rd0 in request.getfixturevalue(corpus):
        i0 = inchi_of(rd0)
        if i0 != inchi_of(to_rdkit(from_rdkit(rd0))):
            bad.append(i0)
    assert bad == []


@mark.parametrize('corpus', ['stereo_sdf', 'nci', 'filter_examples'])
def test_corpus_configurations_match_the_core(corpus, request):
    """The independent anchor: what the core's own InChI writer says about the configuration.

    A `/b` disagreement is allowed only where the core's own layer carries an undefined mark -- one
    record today, a C170 phthalocyanine -- and cannot hide a lost configuration, which produces no `?`.
    """
    bad = []
    for rd0 in request.getfixturevalue(corpus):
        i0, i1 = inchi_of(rd0), molecule_to_inchi(from_rdkit(rd0))
        if layers(i0, 'tms') != layers(i1, 'tms'):
            bad.append(('tms', i0, i1))
        elif layers(i0, 'b') != layers(i1, 'b') and '?' not in layers(i1, 'b').get('b', ''):
            bad.append(('b', i0, i1))
    assert bad == []
