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
"""The structures the version 3 and version 4 fixture corpora are written from, and their answers.

THE BUILDERS ARE THE SPECIFICATION and the committed corpora are the pin: `gen_pach3_corpus.py` writes
each builder's record once, and `test_pach3.py` asserts today's writer still emits those bytes.  Do not
regenerate a corpus to make a test pass -- a differing byte is either a deliberate layout change, in
which case regenerate and say so in the commit, or the defect the fixture exists to catch.

Between them the builders exercise every field of both layouts: elements, an R marker with and without
an index, isotopes, charges either side of zero, a radical, a pinned count, H_UNKNOWN, an aromatic
ring, a dative bond, coordinates, a wedge, all four stereo kinds, an enhanced-stereo group and map
numbers.

Container format: `pach_corpus.load_corpus` reads it, so it is that module's, byte for byte.

PUBLIC COMPOUNDS ONLY: ethanol, benzene, pyridine, aspirin, caffeine, sodium acetate, alanine,
2-butene, an allene, 2-chloro-2'-fluorobiphenyl, butane, dimethyl maleate, muconic acid, BINOL,
fumaric acid.
"""
from pathlib import Path

from chython.core import H_UNKNOWN, MoleculeContainer, STEREO_AND, STEREO_OR, read_smiles
from .pach_corpus import load_corpus


__all__ = ['BUILDERS', 'BGROUP_BUILDERS', 'BGROUP_NAMES', 'V3_PATH', 'V4_PATH', 'BGROUP_PATH',
           'answers', 'answers_bgroup', 'drawn', 'load_corpus']


V3_PATH = Path(__file__).parent / 'pach_v3_corpus.bin.gz'
V4_PATH = Path(__file__).parent / 'pach_v4_corpus.bin.gz'
BGROUP_PATH = Path(__file__).parent / 'pach_bond_group_corpus.bin.gz'

# 2-chloro-2'-fluorobiphenyl.  Built rather than parsed: SMILES has no atropisomer notation.
_BIPHENYL_BONDS = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                   (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1),
                   (0, 6, 1), (1, 12, 1), (7, 13, 1)]
_BIPHENYL_H = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0]

_BUTANE_XY = [(0.0, 0.0), (1.5, 0.0), (2.0, -1.25), (3.5, -1.25)]


def _chlorofluorobiphenyl():
    mol = MoleculeContainer()
    with mol.edit():
        sids = [mol.add_atom(e, implicit_h=h)
                for e, h in zip(['C'] * 12 + ['Cl', 'F'], _BIPHENYL_H)]
        for i, j, o in _BIPHENYL_BONDS:
            mol.add_bond(sids[i], sids[j], o)
    mol.set_parity(sids[0], 1)
    return mol


def _drawn_butane():
    """Four carbons with a drawing, so `has_coordinates` is True and `pack()` chooses version 3.

    Built rather than parsed: a core test may not import `chython.formats`, and `set_xy` registers
    SEG_XY exactly as a molfile read does.
    """
    mol = MoleculeContainer()
    with mol.edit():
        sids = [mol.add_atom('C', implicit_h=3 if i in (0, 3) else 2) for i in range(4)]
        for i in range(3):
            mol.add_bond(sids[i], sids[i + 1], 1)
    for sid, (x, y) in zip(sids, _BUTANE_XY):
        mol.set_xy(sid, x, y)
    return mol


def drawn(mol):
    """A deterministic 2D layout for a molecule that has none, so its version 3 record has a drawing.

    Not a depiction: `clean2d` lives in `interop` and a core test may not import it.  The coordinates
    make no chemical claim, they populate the field version 3 has and version 4 does not, which is the
    whole difference between the two corpora.  A builder that draws its own molecule keeps its drawing.
    """
    if mol.has_coordinates:
        return mol
    for i, k in enumerate(mol.atom_numbers):
        mol.set_xy(k, i * 1.5, (i % 2) * 0.8)
    return mol


def _wedged_alanine():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    ids = list(mol.atom_numbers)
    for i, k in enumerate(ids):
        mol.set_xy(k, i * 1.5, (i % 2) * 0.8)
    mol.set_wedge(ids[1], ids[2], 1)
    return mol


def _grouped_alanine():
    mol = read_smiles('N[C@@H](C)C(=O)O')
    mol.set_stereo_group(mol.atom_numbers[1], STEREO_AND, 3)
    return mol


def _unknown_hydrogens():
    mol = read_smiles('[SeH4]')
    mol.set_hydrogens(mol.atom_numbers[0], H_UNKNOWN)
    return mol


BUILDERS = [
    ('ethanol', lambda: read_smiles('CCO')),
    ('benzene', lambda: read_smiles('c1ccccc1')),
    ('pyridine', lambda: read_smiles('c1ccncc1')),
    ('aspirin', lambda: read_smiles('CC(=O)Oc1ccccc1C(=O)O')),
    ('caffeine', lambda: read_smiles('Cn1cnc2c1c(=O)n(C)c(=O)n2C')),
    ('sodium_acetate', lambda: read_smiles('CC(=O)[O-].[Na+]')),
    ('ammonium_chloride', lambda: read_smiles('[NH4+].[Cl-]')),
    ('ferric_ion', lambda: read_smiles('[Fe+3]')),
    ('isotopes', lambda: read_smiles('[13CH3][18OH]')),
    ('methyl_radical', lambda: read_smiles('[CH3]')),
    ('methanide', lambda: read_smiles('[CH3-]')),
    ('selenium_unknown_h', _unknown_hydrogens),
    ('bare_r_marker', lambda: read_smiles('[R]C')),
    ('indexed_r_marker', lambda: read_smiles('[R7]CC')),
    ('ammine_platinum', lambda: read_smiles('[NH3]->[Pt]')),
    ('alanine_tetrahedral', lambda: read_smiles('N[C@@H](C)C(=O)O')),
    ('alanine_mirror', lambda: read_smiles('N[C@H](C)C(=O)O')),
    ('trans_butene', lambda: read_smiles('C/C=C/C')),
    ('cis_butene', lambda: read_smiles('C/C=C\\C')),
    ('difluorodibromoallene', lambda: read_smiles('FC(Br)=[C@]=C(F)Br')),
    ('chlorofluorobiphenyl', _chlorofluorobiphenyl),
    ('drawn_butane', _drawn_butane),
    ('wedged_alanine', _wedged_alanine),
    ('grouped_alanine', _grouped_alanine),
    ('mapped_methanol', lambda: read_smiles('[CH3:1][OH:2]')),
]


def _maleate_rel1():
    """Dimethyl (2Z)-but-2-enedioate (dimethyl maleate) -- one C=C bond in OR group 1.

    STEBREL1: a single-bond correlation across one double bond.
    """
    mol = read_smiles('COC(=O)/C=C\\C(=O)OC')
    pairs = [(b.n, b.m) for b in mol.bonds()
             if b.order == 2 and mol.atom(b.n).atomic_symbol == 'C'
             and mol.atom(b.m).atomic_symbol == 'C']
    assert len(pairs) == 1, pairs
    with mol.edit() as e:
        e.set_bond_stereo_group(pairs[0][0], pairs[0][1], STEREO_OR, 1)
    return mol


def _muconate_and1():
    """(2E,4E)-hexa-2,4-dienedioic acid -- two C=C bonds in AND group 1.

    STEBRAC1: a two-bond correlation.
    """
    mol = read_smiles('OC(=O)/C=C/C=C/C(=O)O')
    pairs = [(b.n, b.m) for b in mol.bonds()
             if b.order == 2 and mol.atom(b.n).atomic_symbol == 'C'
             and mol.atom(b.m).atomic_symbol == 'C']
    assert len(pairs) == 2, pairs
    with mol.edit() as e:
        for n, m in pairs:
            e.set_bond_stereo_group(n, m, STEREO_AND, 1)
    return mol


def _binol_and1():
    """1,1'-binaphthalene-2,2'-diol -- axial bond in AND group 1.

    STEBRAC1 on the biaryl single bond, the prototypical atropisomer correlation.  The one C-C
    single bond in BINOL is the inter-ring axis; every other C-C bond in the two naphthalene
    systems is aromatic (order 4).
    """
    mol = read_smiles('Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12')
    inter = [(b.n, b.m) for b in mol.bonds()
             if b.order == 1 and mol.atom(b.n).atomic_symbol == 'C'
             and mol.atom(b.m).atomic_symbol == 'C']
    assert len(inter) == 1, inter
    n, m = inter[0]
    with mol.edit() as e:
        e.set_bond_stereo_group(n, m, STEREO_AND, 1)
    return mol


def _alanine_fumarate_one_collection():
    """Alanine's centre and fumarate's axis in ONE AND 1 collection, in a two-component record.

    One namespace in one fixture: a collection is a set of stereocentres, and a centre and an axis are
    both stereocentres, so `stereo_groups()` answers one key with two members and
    `bond_stereo_groups()` answers the axis alone.  Fumarate carries the axis because alanine has no
    stereogenic double bond -- a carbonyl is a bond a group can be STATED on and no unit owns.
    """
    mol = read_smiles('N[C@@H](C)C(=O)O.OC(=O)/C=C/C(=O)O')
    ids = list(mol.atom_numbers)
    mol.set_stereo_group(ids[1], STEREO_AND, 1)
    axis = [p for p, u in mol.chiral_bonds().items() if u['kind'] == 1]
    assert len(axis) == 1, axis
    with mol.edit() as e:
        e.set_bond_stereo_group(axis[0][0], axis[0][1], STEREO_AND, 1)
    return mol


BGROUP_BUILDERS = [
    ('maleate_rel1', _maleate_rel1),
    ('muconate_and1', _muconate_and1),
    ('binol_and1', _binol_and1),
    ('alanine_fumarate_one_collection', _alanine_fumarate_one_collection),
]

BGROUP_NAMES = [name + sfx for name, _ in BGROUP_BUILDERS for sfx in ('_v4', '_v3')]


def answers_bgroup(mol):
    """Enhanced-stereo answers for `BGROUP_BUILDERS`: what the decoder must reproduce.

    Version-independent: an axis's entry is the same in a v3 and a v4 record, so one answer serves
    both of a builder's records.  ONE NAMESPACE, so `stereo_groups()` holds both member spellings and
    `bond_stereo_groups()` is its axis-only view; both are checked, and every member is spelled as a
    list of atom-list positions -- one long for a centre, two for an axis -- so the answers survive a
    renumbering and sort against each other.
    """
    index = {k: i for i, k in enumerate(mol.atom_numbers)}

    def spell(member):
        return sorted(index[n] for n in member) if isinstance(member, tuple) else [index[member]]

    return {
        'groups': {str(k): sorted(spell(m) for m in v) for k, v in mol.stereo_groups().items()},
        'bond_groups': {str(k): sorted(spell(m) for m in v)
                        for k, v in mol.bond_stereo_groups().items()},
    }


def answers(mol, with_drawing):
    """Everything a decode has to reproduce, in ATOM ORDER and by INDEX rather than by stable id.

    The format stores no id -- `pack()` renumbers -- so an answer keyed by one would be asserting
    something the record does not carry.  Indices are what the record's own slots hold.

    `with_drawing` is False for version 4, which has no coordinate block and therefore no wedge either:
    a wedge is a statement about a drawing.  Storing a coordinate the record cannot hold would make the
    fixture disagree with its own decode.
    """
    index = {k: i for i, k in enumerate(mol.atom_numbers)}
    return {
        'atoms': [[a.element, a.r_index, a.isotope, a.charge, int(a.is_radical), a.implicit_h,
                   list(a.xy) if with_drawing and a.xy is not None else None, a.map_number,
                   list(mol.stereo_group_of(a.n))] for a in mol.atoms()],
        'bonds': [[index[b.n], index[b.m], b.order,
                   [index[b.wedge[0]], b.wedge[1]] if with_drawing and b.wedge is not None else None]
                  for b in mol.bonds()],
        'stereo': sorted([u['kind'], index[u['anchor']],
                          sorted(index[r] for r in u['refs'] if r is not None), u['parity']]
                         for u in mol.stereo_units() if u['parity']),
    }
