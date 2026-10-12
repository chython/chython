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
"""Drawings that state no configuration, or state one twice, read from hand-built CTfiles.

Each case is a drawing convention from the IUPAC 2006 graphical-representation recommendations (ST-n
numbers) or the CTfile specification, read through the real parser.  Structures are public compounds:
bromochlorofluoromethane, fluorochlorobromoiodomethane, but-2-ene, 2-chlorobut-2-ene, pent-3-en-2-ol,
3-chlorobutan-2-ol, methyl ethyl sulfoxide and 6,6'-dimethylbiphenyl-2,2'-diol.
"""

import math

from pytest import mark

from ....core import read_smiles
from .._sdf import sniff_version
from .._v2000 import emit_v2000, parse_v2000
from .._v3000 import V3000_STAMP, emit_v3000, parse_v3000


# --- builders -------------------------------------------------------------------------------------

def _v2000(atoms, bonds, *, chiral=0, parities=None, dim=''):
    """`atoms` are ``(symbol, x, y[, z])``, `bonds` ``(a, b, order, stereo)`` 1-based, V2000 codes."""
    parities = parities or {}
    lines = ['probe', f'  {"probe":8}{"0930261200":10}{dim:2}', '',
             f'{len(atoms):3d}{len(bonds):3d}  0  0{chiral:3d}  0  0  0  0  0999 V2000']
    for i, a in enumerate(atoms, 1):
        s, x, y = a[:3]
        z = a[3] if len(a) > 3 else 0.0
        lines.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {s:<3} 0  0{parities.get(i, 0):3d}  0  0  0  0  0  0  0  0  0')
    for a, b, o, st in bonds:
        lines.append(f'{a:3d}{b:3d}{o:3d}{st:3d}')
    lines.append('M  END')
    return lines


def _v3000(atoms, bonds, *, chiral=0, cfg=None, collections=()):
    """As :func:`_v2000`, bond stereo in V3000 ``CFG`` codes: 1 up, 2 either, 3 down."""
    cfg = cfg or {}
    lines = ['probe', '  probe', '', '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB',
             f'M  V30 COUNTS {len(atoms)} {len(bonds)} 0 0 {chiral}', 'M  V30 BEGIN ATOM']
    for i, a in enumerate(atoms, 1):
        s, x, y = a[:3]
        z = a[3] if len(a) > 3 else 0.0
        extra = f' CFG={cfg[i]}' if i in cfg else ''
        lines.append(f'M  V30 {i} {s} {x:.4f} {y:.4f} {z:.4f} 0{extra}')
    lines += ['M  V30 END ATOM', 'M  V30 BEGIN BOND']
    for j, (a, b, o, st) in enumerate(bonds, 1):
        lines.append(f'M  V30 {j} {o} {a} {b}' + (f' CFG={st}' if st else ''))
    lines.append('M  V30 END BOND')
    if collections:
        lines.append('M  V30 BEGIN COLLECTION')
        lines.extend('M  V30 ' + c for c in collections)
        lines.append('M  V30 END COLLECTION')
    lines += ['M  V30 END CTAB', 'M  END']
    return lines


def _read(lines):
    parse = parse_v3000 if sniff_version(lines, []) == V3000_STAMP else parse_v2000
    mol, _, log = parse(lines, []).build()
    return mol, [r.rule for r in log]


def _parity(mol, kind):
    units = [u for u in mol.stereo_units() if u['kind'] == kind]
    assert len(units) == 1, units
    return mol.parity_of(units[0]['anchor'])


def _at(deg, r=1.5):
    return r * math.cos(math.radians(deg)), r * math.sin(math.radians(deg))


def _centre(spec, lens=None):
    """A carbon at the origin with F, Cl, Br[, I] at ``(angle, V2000 stereo)`` from `spec`."""
    lens = lens or [1.5] * len(spec)
    atoms = [('C', 0.0, 0.0)] + [(s, *_at(a, r)) for s, (a, _), r in zip(('F', 'Cl', 'Br', 'I'), spec, lens)]
    bonds = [(1, i + 2, 1, st) for i, (_, st) in enumerate(spec)]
    return atoms, bonds


U, D, W = 1, 6, 4   # V2000 bond stereo: up, down, either

#: (E)-but-2-ene, zigzag; bond 2 is the double bond.
_BUTENE = [('C', 0.0, 0.0), ('C', 0.866, 0.5), ('C', 1.732, 0.0), ('C', 2.598, 0.5)]


def _butene(double=0, first=(1, 2), first_stereo=0, z=None):
    atoms = [a + ((z[i],) if z else ()) for i, a in enumerate(_BUTENE)]
    return atoms, [(*first, 1, first_stereo), (2, 3, 2, double), (3, 4, 1, 0)]


#: A bromochlorofluoromethane tetrahedron in 3D, implicit hydrogen along -z.
_TETRA_3D = [('C', 0.0, 0.0, 0.0), ('F', 0.0, 1.0, 0.4), ('Cl', -0.9, -0.5, 0.4), ('Br', 0.9, -0.5, 0.4)]


# --- 1. either on a double bond, and a wavy bond at its terminal ----------------------------------

def test_a_crossed_double_bond_leaves_the_configuration_unset_v2000():
    """V2000 bond stereo 3 is "cis or trans (either)": the coordinates are not a statement."""
    mol, log = _read(_v2000(*_butene(double=3)))
    assert _parity(mol, 1) == 0
    assert 'wedge:either-bond' in log, log
    mol, log = _read(_v2000(*_butene()))
    assert _parity(mol, 1) == 1, 'the zigzag without the cross is E'


def test_a_crossed_double_bond_leaves_the_configuration_unset_v3000():
    mol, log = _read(_v3000(*_butene(double=2)))
    assert _parity(mol, 1) == 0
    assert 'wedge:either-bond' in log, log


@mark.parametrize('first', [(2, 1), (1, 2)], ids=['narrow-at-terminal', 'narrow-at-methyl'])
def test_a_wavy_bond_at_a_terminal_leaves_the_double_bond_unset(first):
    """IUPAC ST-4.4: a wavy single bond from a double-bond terminal means E/Z unknown."""
    mol, log = _read(_v2000(*_butene(first=first, first_stereo=W)))
    assert _parity(mol, 1) == 0
    assert 'wedge:either-bond' in log, log


def test_a_wavy_bond_owned_by_a_stereocentre_leaves_the_double_bond_alone():
    """(E)-pent-3-en-2-ol with a wavy bond whose narrow end is C2: it is C2's, not the alkene's."""
    atoms = [('C', 0.0, 0.0), ('C', 0.866, 0.5), ('C', 1.732, 0.0), ('C', 2.598, 0.5), ('C', 3.464, 0.0),
             ('O', 0.866, 1.5)]
    bonds = [(1, 2, 1, 0), (2, 3, 1, W), (3, 4, 2, 0), (4, 5, 1, 0), (2, 6, 1, 0)]
    mol, log = _read(_v2000(atoms, bonds))
    assert mol.parity_of(2) == 0
    assert _parity(mol, 1) == 1


def test_either_wins_over_3d_coordinates():
    """In 3D the coordinates state a configuration, and an either bond still says unknown."""
    atoms = list(_TETRA_3D)
    mol, _ = _read(_v2000(atoms, [(1, 2, 1, 0), (1, 3, 1, 0), (1, 4, 1, 0)]))
    assert mol.parity_of(1), 'the tetrahedron alone is configured'
    mol, log = _read(_v2000(atoms, [(1, 2, 1, W), (1, 3, 1, 0), (1, 4, 1, 0)]))
    assert mol.parity_of(1) == 0
    assert 'wedge:either-bond' in log, log

    mol, log = _read(_v2000(*_butene(double=3, z=(0.3, 0.0, 0.0, 0.0))))
    assert _parity(mol, 1) == 0
    assert 'wedge:either-bond' in log, log


# --- 2. either against the atom parity field ------------------------------------------------------

def test_an_either_bond_wins_over_the_parity_field():
    atoms, bonds = _centre([(210, 0), (330, 0), (90, W)])
    mol, log = _read(_v2000(atoms, bonds, parities={1: 1}))
    assert mol.parity_of(1) == 0
    assert 'wedge:either-bond' in log, log
    assert 'wedge:drawing-field-disagree' in log, log

    atoms, bonds = _centre([(210, 0), (330, 0), (90, 0)])
    mol, log = _read(_v2000(atoms, bonds, parities={1: 1}))
    assert mol.parity_of(1), 'a plain drawing leaves the field in charge'


# --- 3. both substituents of one end on the same side ---------------------------------------------

def _chlorobutene(cl):
    atoms = [('C', -0.75, 1.3), ('C', 0.0, 0.0), ('C', 1.5, 0.0), ('C', 2.25, 1.3), ('Cl', *cl)]
    return atoms, [(1, 2, 1, 0), (2, 3, 2, 0), (3, 4, 1, 0), (2, 5, 1, 0)]


def test_two_substituents_on_one_side_leave_the_double_bond_unset():
    """IUPAC ST-4.1: both substituents of a terminal on one side of the double bond is drawn "in no case"."""
    mol, log = _read(_v2000(*_chlorobutene((-1.5, 0.3))))
    assert _parity(mol, 1) == 0
    assert 'wedge:same-side-substituents' in log, log

    mol, log = _read(_v2000(*_chlorobutene((-0.75, -1.3))))
    assert _parity(mol, 1), 'opposite sides are a drawing'
    assert 'wedge:same-side-substituents' not in log, log


def test_the_same_side_verdict_does_not_depend_on_atom_order():
    atoms, bonds = _chlorobutene((-1.5, 0.3))
    order = [5, 4, 3, 2, 1]                      # new position -> old atom number
    new = {old: i + 1 for i, old in enumerate(order)}
    atoms = [atoms[old - 1] for old in order]
    bonds = [(new[a], new[b], o, s) for a, b, o, s in bonds]
    mol, log = _read(_v2000(atoms, bonds))
    assert _parity(mol, 1) == 0
    assert 'wedge:same-side-substituents' in log, log


# --- 4. near-degenerate tetrahedral geometry ------------------------------------------------------

#: IUPAC ST-1 drawings that state a configuration, and the parity each reads as.
_ACCEPTABLE = {
    '1.1.1': ([(210, 0), (330, 0), (60, U), (120, D)], None, 2),
    '1.1.2': ([(90, 0), (210, 0), (330, 0), (30, U)], None, 1),
    '1.1.3': ([(0, U), (90, D), (180, U), (270, 0)], None, 2),
    '1.1.4': ([(0, 0), (90, 0), (180, 0), (270, U)], None, 1),
    '1.1.5': ([(90, U), (210, U), (330, U), (30, 0)], None, 2),
    '1.2.1': ([(210, 0), (330, 0), (90, U)], None, 2),
    '1.2.1 hashed': ([(210, 0), (330, 0), (90, D)], None, 1),
    '1.2.2': ([(90, 0), (0, U), (300, D)], None, 2),
    '1.2.3': ([(210, 0), (330, 0), (270, U)], None, 2),
    '1.2.3 short wedge': ([(210, 0), (330, 0), (270, U)], [1.5, 1.5, 4.5], 2),
    '1.2.3 long plain': ([(210, 0), (330, 0), (270, U)], [4.5, 4.5, 1.0], 2),
    '1.2.4': ([(0, 0), (60, 0), (120, U)], None, 2),
    '1.2.4 short wedge': ([(0, 0), (60, 0), (120, U)], [3, 3, 1], 2),
    '1.2.5': ([(90, 0), (210, U), (330, U)], None, 2),
}


@mark.parametrize('name', list(_ACCEPTABLE))
def test_an_acceptable_drawing_reads_its_configuration(name):
    spec, lens, parity = _ACCEPTABLE[name]
    mol, log = _read(_v2000(*_centre(spec, lens)))
    assert mol.parity_of(1) == parity, log
    assert 'wedge:ambiguous-drawing' not in log, log


#: Drawings a few degrees off a degenerate layout: the sign is decided by those degrees.
_NEAR_DEGENERATE = {
    'square planar +1': [(0, 0), (181, 0), (90, U), (270, D)],
    'two up two down +1': [(45, U), (136, U), (225, D), (315, D)],
    'T-shape 176': [(0, 0), (176, 0), (90, U)],
    'T-shape 179': [(0, 0), (179, 0), (90, U)],
    'T-shape 181': [(0, 0), (181, 0), (90, U)],
    'trigonal planar +1': [(91, 0), (210, U), (330, D)],
}


@mark.parametrize('name', list(_NEAR_DEGENERATE))
def test_a_near_degenerate_drawing_is_left_unset(name):
    mol, log = _read(_v2000(*_centre(_NEAR_DEGENERATE[name])))
    assert mol.parity_of(1) == 0
    assert 'wedge:ambiguous-drawing' in log, log


# --- 5. drawings dropped without a word -----------------------------------------------------------

def test_a_wedge_on_a_double_bond_is_reported():
    mol, log = _read(_v2000(*_butene(double=U)))
    assert 'wedge:on-double-bond' in log, log
    assert _parity(mol, 1) == 1, 'the double bond itself still reads from the coordinates'


def test_a_sulfoxide_wedge_on_s_o_is_read_and_not_reported():
    atoms = [('C', 0.0, 0.0), ('S', 0.866, 0.5), ('O', 0.866, 2.0), ('C', 1.732, 0.0), ('C', 2.598, 0.5)]
    bonds = [(1, 2, 1, 0), (2, 3, 2, U), (2, 4, 1, 0), (4, 5, 1, 0)]
    mol, log = _read(_v2000(atoms, bonds))
    assert mol.parity_of(2)
    assert 'wedge:on-double-bond' not in log, log


#: 3-chlorobutan-2-ol, zigzag: C2 carries O5, C3 carries Cl6.
_CHLOROBUTANOL = [('C', 0.0, 0.0), ('C', 0.866, 0.5), ('C', 1.732, 0.0), ('C', 2.598, 0.5),
                  ('O', 0.866, 1.5), ('Cl', 1.732, -1.0)]


@mark.parametrize('code, reported', [(D, True), (U, False)], ids=['hashed', 'solid'])
def test_a_hashed_bond_between_two_centres_is_reported(code, reported):
    """IUPAC ST-0.5: a hashed bond between two stereocentres is ambiguous; it is read at the narrow end."""
    bonds = [(1, 2, 1, 0), (2, 3, 1, code), (3, 4, 1, 0), (2, 5, 1, 0), (3, 6, 1, 0)]
    mol, log = _read(_v2000(_CHLOROBUTANOL, bonds))
    assert mol.parity_of(2)
    assert mol.parity_of(3) == 0
    assert ('wedge:hash-between-centres' in log) is reported, log


def test_a_wedge_contradicting_3d_coordinates_is_reported():
    plain, _ = _read(_v2000(_TETRA_3D, [(1, 2, 1, 0), (1, 3, 1, 0), (1, 4, 1, 0)]))
    parity = plain.parity_of(1)
    reported = []
    for code in (U, D):
        mol, log = _read(_v2000(_TETRA_3D, [(1, 2, 1, code), (1, 3, 1, 0), (1, 4, 1, 0)]))
        assert mol.parity_of(1) == parity, 'the coordinates win'
        reported.append('wedge:3d-wedge-overridden' in log)
    assert sorted(reported) == [False, True], 'only the contradicting wedge is reported'


@mark.parametrize('dim, z, reported', [('3D', 0.0, True), ('2D', 0.4, True), ('2D', 0.0, False),
                                       ('3D', 0.4, False), ('', 0.4, False)])
def test_a_dimension_stamp_disagreeing_with_z_is_reported(dim, z, reported):
    atoms = [('C', 0.0, 0.0, z), ('O', 1.5, 0.0, 0.0)]
    _, log = _read(_v2000(atoms, [(1, 2, 1, 0)], dim=dim))
    assert ('ctab:dimensionality-mismatch' in log) is reported, log


#: 6,6'-dimethylbiphenyl-2,2'-diol: pivots 8 and 9; C2 (ortho to 8) carries the methyl C1.
_BIARYL = [(-0.3572, -1.4438), (0.3572, -1.0312), (1.0717, -1.4437), (1.7862, -1.0312), (1.7862, -0.2062),
           (1.0717, 0.2063), (1.0717, 1.0313), (0.3572, -0.2062), (-0.3572, 0.2062), (-1.0717, -0.2063),
           (-1.0717, -1.0313), (-1.7862, 0.2062), (-1.7862, 1.0312), (-1.0717, 1.4437), (-0.3572, 1.0312),
           (0.3572, 1.4438)]
_BIARYL_SYMBOLS = 'CCCCCCOCCCCCCCCO'
_BIARYL_BONDS = [(1, 2, 1), (2, 3, 4), (8, 2, 4), (3, 4, 4), (4, 5, 4), (5, 6, 4), (6, 7, 1), (6, 8, 4),
                 (8, 9, 1), (9, 10, 4), (9, 15, 4), (10, 11, 1), (10, 12, 4), (12, 13, 4), (13, 14, 4),
                 (14, 15, 4), (15, 16, 1)]


def _biaryl(wedges):
    atoms = [(s, x, y) for s, (x, y) in zip(_BIARYL_SYMBOLS, _BIARYL)]
    bonds = []
    for a, b, o in _BIARYL_BONDS:
        if (b, a) in wedges:
            a, b = b, a
        bonds.append((a, b, o, wedges.get((a, b), 0)))
    return _v2000(atoms, bonds)


def test_a_ring_bond_wedge_configures_the_axis():
    mol, log = _read(_biaryl({(8, 2): U}))
    assert _parity(mol, 3)
    assert 'wedge:atropisomer-unsupported-mark' not in log, log


@mark.parametrize('wedge', [(2, 1), (8, 9)], ids=['ortho-substituent', 'pivot-bond'])
def test_an_atropisomer_mark_off_the_ring_bonds_is_reported(wedge):
    """IUPAC ST-1.7: the axis is marked on a ring bond at a pivot, not on a substituent or the pivot bond."""
    mol, log = _read(_biaryl({wedge: U}))
    assert _parity(mol, 3) == 0
    assert 'wedge:atropisomer-unsupported-mark' in log, log


# --- 7. V3000 collections -------------------------------------------------------------------------

_BONDS_CB = [(1, 2, 1, 0), (2, 3, 1, 0), (3, 4, 1, 0), (2, 5, 1, 1), (3, 6, 1, 0)]


def test_an_unconfigured_collection_member_is_reported():
    mol, log = _read(_v3000(_CHLOROBUTANOL, _BONDS_CB, collections=['MDLV30/STERAC1 ATOMS=(2 2 3)']))
    assert mol.parity_of(3) == 0
    assert 'v3000:collection-member-unconfigured' in log, log
    assert mol.stereo_group_of(2) == (3, 1)


def test_an_atom_in_two_collections_keeps_the_first():
    lines = _v3000(_CHLOROBUTANOL, _BONDS_CB,
                   collections=['MDLV30/STEREL1 ATOMS=(1 2)', 'MDLV30/STERAC2 ATOMS=(1 2)'])
    mol, log = _read(lines)
    assert 'v3000:collection-duplicate-member' in log, log
    assert mol.stereo_group_of(2) == (2, 1), 'the first collection naming the atom wins'


# --- 9. the chiral flag ---------------------------------------------------------------------------

@mark.parametrize('version', ['v2000', 'v3000'])
@mark.parametrize('chiral', [0, 1])
def test_the_chiral_flag_1_is_reported_ignored(version, chiral):
    build = _v2000 if version == 'v2000' else _v3000
    mol, log = _read(build(_CHLOROBUTANOL, _BONDS_CB, chiral=chiral))
    assert mol.parity_of(2)
    assert mol.stereo_group_of(2) == (0, 0), 'the flag is not an ABS label'
    assert ('ctab:chiral-flag-ignored' in log) is bool(chiral), log


def test_the_chiral_flag_without_a_configuration_is_silent():
    bonds = [(a, b, o, 0) for a, b, o, _ in _BONDS_CB]
    _, log = _read(_v2000(_CHLOROBUTANOL, bonds, chiral=1))
    assert 'ctab:chiral-flag-ignored' not in log, log


@mark.parametrize('tag', ['STERAC1', 'STEREL1'])
def test_a_chiral_flag_contradicting_a_collection_is_reported(tag):
    lines = _v3000(_CHLOROBUTANOL, _BONDS_CB, chiral=1, collections=[f'MDLV30/{tag} ATOMS=(1 2)'])
    _, log = _read(lines)
    assert 'v3000:chiral-flag-contradicts-collection' in log, log
    lines = _v3000(_CHLOROBUTANOL, _BONDS_CB, chiral=1, collections=['MDLV30/STEABS ATOMS=(1 2)'])
    _, log = _read(lines)
    assert 'v3000:chiral-flag-contradicts-collection' not in log, log


@mark.parametrize('emit, rule', [(emit_v2000, 'v2000:unlabelled-sign-flag'),
                                 (emit_v3000, 'v3000:unlabelled-sign-flag')], ids=['v2000', 'v3000'])
def test_writing_the_flag_for_unlabelled_signs_says_what_readers_assume(emit, rule):
    mol = read_smiles('C[C@H](O)CC')
    log = []
    lines, _ = emit(mol, log=log)
    counts = lines[3] if emit is emit_v2000 else next(x for x in lines if 'COUNTS' in x)
    assert counts.split()[-1 if emit is emit_v3000 else 4] == '1', counts
    assert rule in [r.rule for r in log], log

    labelled = read_smiles('C[C@H](O)CC |a:1|')
    log = []
    emit_v3000(labelled, log=log)
    assert 'v3000:unlabelled-sign-flag' not in [r.rule for r in log], log


# --- 10. atropisomer axes in collections, on write ------------------------------------------------

def _grouped_biaryl():
    mol, _ = _read(_biaryl({(8, 2): U}))
    (pair,) = list(mol.chiral_bonds())
    with mol.edit():
        mol.set_bond_stereo_group(pair[0], pair[1], 3, 1)
    return mol, pair


def test_an_axis_collection_is_written_on_its_owner_atom():
    mol, pair = _grouped_biaryl()
    anchor = next(u['anchor'] for u in mol.stereo_units() if u['kind'] == 3)
    lines, _ = emit_v3000(mol)
    position = list(mol.atom_numbers).index(anchor) + 1
    assert any(f'MDLV30/STERAC1 ATOMS=(1 {position})' in x for x in lines), lines
    assert not any('STEB' in x for x in lines), lines
    back, _ = _read(lines)
    assert back.bond_stereo_groups() == {(3, 1): [pair]}, back.stereo_groups()


def test_a_double_bond_collection_keeps_the_bond_spelling():
    mol = read_smiles('C[C@H](O)/C=C/C')
    (pair,) = list(mol.chiral_bonds())
    with mol.edit():
        mol.set_bond_stereo_group(pair[0], pair[1], 3, 1)
    lines, _ = emit_v3000(mol)
    assert any('MDLV30/STEBRAC1 BONDS=' in x for x in lines), lines
