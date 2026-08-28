# -*- coding: utf-8 -*-
from chython import smiles, MoleculeContainer


def test_smiles_parse_ethanol():
    mol = smiles('CCO')
    assert len(mol) == 3
    assert str(mol) == 'CCO'


def test_smiles_parse_benzene():
    mol = smiles('c1ccccc1')
    assert len(mol) == 6
    assert str(mol) == 'c1ccccc1'


def test_molecule_equality():
    mol1 = smiles('CCO')
    mol2 = smiles('OCC')
    assert mol1 == mol2
    assert hash(mol1) == hash(mol2)


def test_molecule_inequality():
    mol1 = smiles('CCO')
    mol2 = smiles('CC=O')
    assert mol1 != mol2


def test_brutto():
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
    assert mol.brutto == {'C': 9, 'H': 8, 'O': 4}


def test_molecular_mass():
    mol = smiles('O')  # water
    mass = mol.molecular_mass
    assert 16 < mass < 20


def test_rings():
    mol = smiles('c1ccc2ccccc2c1')  # naphthalene
    assert mol.rings_count == 2
    assert len(mol.sssr) == 2


def test_aromatic_rings():
    mol = smiles('c1ccccc1')
    assert len(mol.aromatic_rings) == 1


def test_connected_components():
    mol = smiles('[Cl-].[Na+]')
    assert mol.connected_components_count == 2
    parts = mol.split()
    assert len(parts) == 2


def test_add_atom_bond():
    mol = MoleculeContainer()
    n1 = mol.add_atom('C')
    n2 = mol.add_atom('C')
    n3 = mol.add_atom('O')
    mol.add_bond(n1, n2, 1)
    mol.add_bond(n2, n3, 1)
    assert str(mol) == 'CCO'


def test_implicit_hydrogens():
    mol = smiles('C')  # methane
    atom = mol.atom(next(iter(mol)))
    assert atom.implicit_hydrogens == 4


def test_copy():
    mol = smiles('CCO')
    mol2 = mol.copy()
    assert str(mol) == str(mol2)
    assert mol is not mol2


def test_substructure():
    mol = smiles('Cc1ccccc1')  # toluene
    atoms = list(mol)
    ring_atoms = atoms[1:]
    sub = mol.substructure(ring_atoms)
    assert len(sub) == 6


def test_union():
    mol1 = smiles('C')
    mol2 = smiles('O')
    merged = mol1.union(mol2, remap=True)
    assert len(merged) == 2


def test_explicify_implicify():
    mol = smiles('C')
    n_before = len(mol)
    mol.explicify_hydrogens()
    assert len(mol) == 5  # C + 4H
    mol.implicify_hydrogens()
    assert len(mol) == n_before


def test_pickle():
    from pickle import loads, dumps
    mol = smiles('c1ccccc1')
    data = dumps(mol)
    mol2 = loads(data)
    assert str(mol) == str(mol2)


def test_format_specifiers():
    mol = smiles('[CH3:1][OH:2]')
    s_map = format(mol, 'm')
    assert ':1]' in s_map or ':2]' in s_map


def test_split_metal_salts_refreshes_labels():
    """Breaking the metal bond must leave the acceptor's marks matching its real environment.

    flush_cache() clears cached properties only, so without calc_labels() the anion keeps
    counting the metal that just left: it reports one neighbor too many and one heteroatom too
    many, and every D/x query on it fails.
    """
    mol = smiles('c1ccccc1C(=O)O[K]')  # potassium benzoate, metal written as a covalent bond
    assert mol.split_metal_salts()

    for n, atom in mol.atoms():
        assert atom.neighbors == sum(b != 8 for b in mol._bonds[n].values()), f'atom {n}'

    anion = next(n for n, a in mol.atoms() if a.charge == -1)
    assert mol.atom(anion).neighbors == 1
    assert mol.atom(anion).heteroatoms == 0


def test_split_metal_salts_then_neutralize():
    """The stale marks made neutralize() a no-op: its base rules never saw the carboxylate."""
    covalent = smiles('c1ccccc1C(=O)O[K]')
    covalent.split_metal_salts()
    ionic = smiles('c1ccccc1C(=O)[O-].[K+]')

    assert covalent == ionic
    assert covalent.neutralize(keep_charge=False)
    assert ionic.neutralize(keep_charge=False)
    assert covalent == ionic
    # the proton is back on the acid; the cation is remove_metals()' business, not neutralize()'s
    assert not any(a.charge < 0 for _, a in covalent.atoms())
