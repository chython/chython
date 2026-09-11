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
from itertools import permutations
from random import Random

import pytest
from chython.core import AutomorphismBudgetExceeded, MoleculeContainer


_SYMBOL = re.compile(r'[A-Z][a-z]?')


def _symbols(atoms):
    """'CFClBr' -> ['C', 'F', 'Cl', 'Br'].

    An element symbol is one capital optionally followed by one lowercase, so this split is
    exact. Do NOT iterate the string directly: 'Cl' is two characters.
    """
    return _SYMBOL.findall(atoms)


def _mol(*, atoms, bonds):
    """atoms: concatenated element symbols. bonds: (i, j, order) over 0-based positions."""
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e) for e in _symbols(atoms)]
        for i, j, o in bonds:
            m.add_bond(sids[i], sids[j], o)
    return m, sids


def _edges(m):
    """Every bond once, as (sid_a, sid_b, order) with sid_a < sid_b."""
    return {(min(a, b), max(a, b), m.order_of(a, b))
            for a in m.atom_numbers for b in m.neighbors_of(a)}


def test_benzene_orbit_is_one_class():
    # Kekule benzene, chosen deliberately now that order 4 IS storable: the aromatic spelling makes
    # all six bonds alike and is the easy case, while a fixed Kekule ring is LESS symmetric and is
    # the one that can go wrong. The alternating C6 still has the rotation by two and the reflections
    # that keep all six carbons in one orbit, so the orbit count is the same for a different reason.
    m, sids = _mol(atoms='C' * 6,
                   bonds=[(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)])
    orbits = m.automorphism_orbits()
    assert len(set(orbits.values())) == 1
    assert not m.is_asymmetric()


def test_ethanol_is_asymmetric():
    m, sids = _mol(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])
    assert m.is_asymmetric()
    assert len(set(m.automorphism_orbits().values())) == 3


def test_isobutane_methyls_share_an_orbit():
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    orbits = m.automorphism_orbits()
    assert orbits[sids[1]] == orbits[sids[2]] == orbits[sids[3]]
    assert orbits[sids[0]] != orbits[sids[1]]


def test_bond_order_is_preserved_by_automorphisms():
    # 1,3-butadiene: C=C-C=C. The two ends swap; the C=C and C-C bonds must not be confused.
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 2), (1, 2, 1), (2, 3, 2)])
    orbits = m.automorphism_orbits()
    assert orbits[sids[0]] == orbits[sids[3]]
    assert orbits[sids[1]] == orbits[sids[2]]
    assert orbits[sids[0]] != orbits[sids[1]]


def test_seeded_colouring_breaks_symmetry():
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    # Marking one methyl by hand takes the group from S3 down to the swap of the other two.
    seed = {sid: (1 if sid != sids[1] else 2) for sid in sids}
    assert not m.is_asymmetric(seed=seed)
    # Marking all three collapses it to the identity: a seed reaches the group, not just the
    # partition.
    seed = {sid: pos for pos, sid in enumerate(sids)}
    assert m.is_asymmetric(seed=seed)
    assert not m.is_asymmetric()


def _brute_force_orbits(m):
    """Orbits by enumerating every permutation of the stable ids -- the oracle, O(n!)."""
    sids = list(m.atom_numbers)
    edges = _edges(m)
    # The whole atom colour, not just the element: an oracle that compares less than
    # `_atom_colour_equal` does cannot catch a regression in the fields it leaves out.
    colour = {s: (m.element_of(s), m.charge_of(s), m.isotope_of(s), m.radical_of(s),
                  m.implicit_h_of(s)) for s in sids}
    parent = {s: s for s in sids}

    def find(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    for perm in permutations(sids):
        sigma = dict(zip(sids, perm))
        if any(colour[s] != colour[t] for s, t in sigma.items()):
            continue
        if {(min(sigma[a], sigma[b]), max(sigma[a], sigma[b]), o) for a, b, o in edges} != edges:
            continue
        for s, t in sigma.items():
            ra, rb = find(s), find(t)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for s in sids:
        groups.setdefault(find(s), set()).add(s)
    return {frozenset(g) for g in groups.values()}


def _core_orbits(m):
    groups = {}
    for sid, orbit in m.automorphism_orbits().items():
        groups.setdefault(orbit, set()).add(sid)
    return {frozenset(g) for g in groups.values()}


def test_orbits_agree_with_brute_force():
    # The search verifies candidates exactly, so it must agree with plain enumeration. Small
    # public structures only: the oracle is factorial.
    cases = [dict(atoms='C' * 5, bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)]),  # pentane
             dict(atoms='C' * 6,                                                       # cyclohexane
                  bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)]),
             dict(atoms='CClClClCl',                                                   # CCl4
                  bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)]),
             dict(atoms='COCOO', bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 2), (2, 4, 1)]),  # MeO-CO-OH
             dict(atoms='CCCCCC', bonds=[(0, 1, 1), (2, 3, 1), (4, 5, 1)]),            # 3 ethanes
             dict(atoms='NCCNCC',                                                      # 2 x ethylamine
                  bonds=[(0, 1, 1), (1, 2, 1), (3, 4, 1), (4, 5, 1)])]
    for case in cases:
        m, _ = _mol(**case)
        assert _core_orbits(m) == _brute_force_orbits(m), case['atoms']

    # And the same against structures whose symmetry turns on the atom-record fields the plain
    # `_mol` builder cannot set: 2-13C-propane (isotope), the glycine zwitterion (charge), and
    # 2-propyl radical (radical plus hydrogen count). Each is the symmetric skeleton with one
    # field broken, so the oracle only agrees if the search reads that field too.
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=3), m.add_atom('C', implicit_h=2, isotope=13),
                m.add_atom('C', implicit_h=3)]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 1)
    assert _core_orbits(m) == _brute_force_orbits(m)

    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('N', implicit_h=3, charge=1), m.add_atom('C', implicit_h=2),
                m.add_atom('C'), m.add_atom('O', charge=-1), m.add_atom('O')]
        for i, j, o in [(0, 1, 1), (1, 2, 1), (2, 3, 1), (2, 4, 2)]:
            m.add_bond(sids[i], sids[j], o)
    assert _core_orbits(m) == _brute_force_orbits(m)

    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=3), m.add_atom('C', implicit_h=1, radical=True),
                m.add_atom('C', implicit_h=3)]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 1)
    assert _core_orbits(m) == _brute_force_orbits(m)


def test_orbits_are_finer_than_symmetry_ranks():
    # Cyclopropane and cyclobutane in one record. Refinement cannot split them -- every atom in
    # either ring is a carbon with two ring neighbours of its own class -- so atoms_order reports
    # one class, while no automorphism can map a three-ring atom onto a four-ring one. Orbits are
    # a subdivision of the ranks, and this is the case that proves the verification does the work.
    m, sids = _mol(atoms='C' * 7,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 0, 1),
                          (3, 4, 1), (4, 5, 1), (5, 6, 1), (6, 3, 1)])
    assert m.atoms_order_classes == 1
    orbits = m.automorphism_orbits()
    assert len({orbits[s] for s in sids[:3]}) == 1
    assert len({orbits[s] for s in sids[3:]}) == 1
    assert orbits[sids[0]] != orbits[sids[3]]


def test_many_identical_fragments_still_share_orbits():
    # Six water molecules in one record: two orbits, the oxygens and the hydrogens. The group has
    # 6! * 2**6 members, far more than any row budget would hold, which is why the search asks
    # per unresolved pair instead of enumerating the group -- an enumeration truncated by a cap
    # returns the permutations of the last few atoms only and reports five orbits here.
    bonds = []
    for k in range(6):
        bonds += [(3 * k, 3 * k + 1, 1), (3 * k, 3 * k + 2, 1)]
    m, sids = _mol(atoms='OHH' * 6, bonds=bonds)
    orbits = m.automorphism_orbits()
    assert len(set(orbits.values())) == 2
    assert len({orbits[s] for s in sids[::3]}) == 1


def _cycloalkanes(sizes):
    """One record holding a cycloalkane per entry of `sizes`, in the order given.

    All-carbon rings with two implicit hydrogens each, i.e. cyclopropane, cyclobutane and friends.
    """
    m = MoleculeContainer()
    with m.edit():
        for k in sizes:
            ring = [m.add_atom('C', implicit_h=2) for _ in range(k)]
            for i in range(k):
                m.add_bond(ring[i], ring[(i + 1) % k], 1)
    return m


def test_orbits_do_not_depend_on_the_order_fragments_were_built_in():
    # 11 cyclopropanes and 11 cyclobutanes in one record: 77 atoms, two orbits -- every CH2 of a
    # three-ring is equivalent to every other, likewise the four-rings, and no automorphism
    # crosses between them. Ordinary chemistry, and well inside any budget.
    #
    # This is the regression guard for the budget being PER SEARCH. A single budget shared by
    # every pair search in the call ran out here and reported 45 orbits for the blocked order --
    # four separate orbits for one cyclobutane's four carbons, which is four invented
    # stereocentres -- while the interleaved order stayed at 2. Orbits must be a function of the
    # structure, so both orders must agree, and they must agree on the truth.
    for sizes in ([3] * 11 + [4] * 11, [3, 4] * 11):
        m = _cycloalkanes(sizes)
        orbits = m.automorphism_orbits()
        assert len(orbits) == 77, sizes
        assert len(set(orbits.values())) == 2, sizes


def test_orbits_survive_many_more_fragments_than_that():
    # The same shape at far more fragments, plus a third ring size so the answer is not trivially
    # two: 20 x cyclopropane + 20 x cyclobutane + 20 x cyclopentane.
    #
    # No test asserts the truncation path itself (AutomorphismBudgetExceeded from
    # automorphism_orbits, False from is_asymmetric), because nothing that fits in a test reaches
    # it: with candidates drawn from adjacency, this record at 2800 atoms still returns exact
    # orbits inside the per-search budget, and the budgets are what a pathological graph -- not a
    # compound -- would need. Both paths were checked by building with the two constants lowered.
    m = _cycloalkanes([3, 4, 5] * 20)
    assert len(set(m.automorphism_orbits().values())) == 3
    m = _cycloalkanes([3] * 20 + [4] * 20 + [5] * 20)
    assert len(set(m.automorphism_orbits().values())) == 3


def _flat_seed_orbits(m):
    """Orbits under a seed that says every atom is alike.

    That empties the refinement of everything but topology -- round 0 is one class, and later
    rounds only fold in neighbour classes and bond orders -- so `cls` no longer knows an element
    from an element and the atom colour comparison inside the search is the ONLY thing left that
    can separate two atoms of equal topology.
    """
    return m.automorphism_orbits(seed=dict.fromkeys(m.atom_numbers, 1))


def test_atom_colour_is_checked_when_the_seed_cannot_tell_atoms_apart():
    # Every case is a symmetric skeleton whose two ends differ in exactly one atom-record field.
    # Under a flat seed the refinement puts those two ends in one class, so an orbit count of 3
    # (or 4) means the search's own colour comparison rejected the swap. Delete a field from
    # `_atom_colour_equal` and one of these drops to 2 (or 3).
    m, sids = _mol(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])        # ethanol: element
    assert len(set(_flat_seed_orbits(m).values())) == 3

    m, sids = _mol(atoms='FCF', bonds=[(0, 1, 1), (1, 2, 1)])        # difluoromethane: control
    assert len(set(_flat_seed_orbits(m).values())) == 2

    m = MoleculeContainer()                                          # 1-13C-propane: isotope
    with m.edit():
        sids = [m.add_atom('C', implicit_h=3, isotope=13), m.add_atom('C', implicit_h=2),
                m.add_atom('C', implicit_h=3)]
        m.add_bond(sids[0], sids[1], 1)
        m.add_bond(sids[1], sids[2], 1)
    assert len(set(_flat_seed_orbits(m).values())) == 3

    m = MoleculeContainer()                        # propane-1,3-diyl anion/cation: charge, and
    with m.edit():                                 # charge ALONE -- both ends carry two hydrogens,
        sids = [m.add_atom('C', implicit_h=2, charge=-1),   # so every other field is equal and the
                m.add_atom('C', implicit_h=2),             # charge comparison is the only thing
                m.add_atom('C', implicit_h=2, charge=1)]    # that can refuse to swap them. A
        for i in range(2):                                  # minimal witness for the comparator,
            m.add_bond(sids[i], sids[i + 1], 1)             # not a compound anyone would isolate.
    assert len(set(_flat_seed_orbits(m).values())) == 3

    m = MoleculeContainer()                                          # propane vs its 1-radical
    with m.edit():                                                   # : radical, hydrogen count
        sids = [m.add_atom('C', implicit_h=2, radical=True), m.add_atom('C', implicit_h=2),
                m.add_atom('C', implicit_h=2)]
        for i in range(2):
            m.add_bond(sids[i], sids[i + 1], 1)
    assert len(set(_flat_seed_orbits(m).values())) == 3
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=3), m.add_atom('C', implicit_h=2),
                m.add_atom('C', implicit_h=2)]
        for i in range(2):
            m.add_bond(sids[i], sids[i + 1], 1)
    assert len(set(_flat_seed_orbits(m).values())) == 3

    # And a coarse seed must not lose symmetry that is really there: isobutane's three methyls
    # stay one orbit, so the colour comparison is rejecting only what it should.
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    orbits = _flat_seed_orbits(m)
    assert len(set(orbits.values())) == 2
    assert orbits[sids[1]] == orbits[sids[2]] == orbits[sids[3]]


def test_empty_and_single_atom_are_asymmetric():
    m = MoleculeContainer()
    assert m.is_asymmetric()
    assert m.automorphism_orbits() == {}
    m, sids = _mol(atoms='C', bonds=[])
    assert m.is_asymmetric()
    assert m.automorphism_orbits() == {sids[0]: 1}


def test_orbits_are_one_based_and_dense():
    m, sids = _mol(atoms='C' * 5, bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)])
    values = sorted(set(m.automorphism_orbits().values()))
    assert values == list(range(1, len(values) + 1))


def test_orbits_are_keyed_by_stable_id():
    m, sids = _mol(atoms='OCO', bonds=[(0, 1, 1), (1, 2, 1)])
    assert set(m.automorphism_orbits()) == set(sids)
    with m.edit():
        m.delete_atom(sids[0])
    assert set(m.automorphism_orbits()) == {sids[1], sids[2]}


def test_seed_is_read_by_orbits_too():
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    seed = {sid: (1 if sid != sids[1] else 2) for sid in sids}
    orbits = m.automorphism_orbits(seed=seed)
    assert len(set(orbits.values())) == 3
    assert orbits[sids[2]] == orbits[sids[3]]  # the two methyls left alike stay equivalent
    assert orbits[sids[1]] != orbits[sids[2]]


def test_seed_conversion_is_the_one_refined_order_uses():
    # the same rejections, because the same conversion runs
    m, sids = _mol(atoms='CCC', bonds=[(0, 1, 1), (1, 2, 1)])
    with pytest.raises(KeyError):
        m.automorphism_orbits(seed={sids[0]: 1, sids[1]: 1})
    with pytest.raises(KeyError):
        m.is_asymmetric(seed={sids[0]: 1, sids[1]: 1})
    with pytest.raises(OverflowError):
        m.automorphism_orbits(seed=dict.fromkeys(sids, -1))
    with pytest.raises(OverflowError):
        m.is_asymmetric(seed=dict.fromkeys(sids, -1))


def _relabel(m, rng):
    """Rebuild m with its atoms added in a random order: same molecule, different slots."""
    perm = list(m.atom_numbers)
    rng.shuffle(perm)
    out = MoleculeContainer()
    with out.edit():
        fresh = {sid: out.add_atom(m.element_of(sid)) for sid in perm}
        for a, b, order in _edges(m):
            out.add_bond(fresh[a], fresh[b], order)
    return out


def _canonical_string(m):
    """A labelling-independent encoding: atoms and bonds keyed by canonical position."""
    pos = m.canonical_order()
    atoms = [None] * len(pos)
    for sid, p in pos.items():
        atoms[p] = m.element_of(sid)
    edges = sorted((min(pos[a], pos[b]), max(pos[a], pos[b]), order)
                   for a, b, order in _edges(m))
    return repr((atoms, edges))


def test_canonical_order_is_relabeling_invariant():
    # cubane skeleton: 8 equivalent carbons, the symmetry a relabeling-variant order oscillates on
    m, sids = _mol(atoms='C' * 8,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
                          (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
                          (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)])
    rng = Random(20260901)
    strings = {_canonical_string(_relabel(m, rng)) for _ in range(60)}
    assert len(strings) == 1, f'{len(strings)} distinct canonical forms'


def test_canonical_order_is_a_permutation():
    m, sids = _mol(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])
    assert sorted(m.canonical_order().values()) == [0, 1, 2]


def test_asymmetric_molecule_needs_no_search():
    m, sids = _mol(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])
    order = m.canonical_order()
    assert len(set(order.values())) == 3


# --- beyond the brief: the same invariance on the shapes that stress it differently -------------

_INVARIANCE_CASES = [
    # adamantane: 10 carbons, 3 orbits, a cage the refinement alone cannot make discrete
    ('adamantane', dict(atoms='C' * 10,
                        bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                               (0, 6, 1), (2, 7, 1), (4, 8, 1), (6, 9, 1), (7, 9, 1), (8, 9, 1)])),
    # Kekule benzene: symmetry that survives a fixed alternation, so bond order has to be in
    # the certificate for the two Kekule-inequivalent positions to stay apart
    ('benzene', dict(atoms='C' * 6,
                     bonds=[(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1)])),
    # three ethanes: disconnected, so the branch order is over components as well as atoms
    ('three ethanes', dict(atoms='C' * 6, bonds=[(0, 1, 1), (2, 3, 1), (4, 5, 1)])),
    # cyclopropane + cyclobutane: one refinement class, two orbits -- the case where the
    # partition cannot help and the search has to
    ('C3 + C4 rings', dict(atoms='C' * 7,
                           bonds=[(0, 1, 1), (1, 2, 1), (2, 0, 1),
                                  (3, 4, 1), (4, 5, 1), (5, 6, 1), (6, 3, 1)])),
    # naphthalene skeleton: fused rings, 3 orbits
    ('naphthalene', dict(atoms='C' * 10,
                         bonds=[(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1),
                                (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 0, 1), (4, 9, 1)])),
    # six waters: 6! * 2**6 automorphisms, the record that broke a row-capped group enumeration
    ('six waters', dict(atoms='OHH' * 6,
                        bonds=[b for k in range(6)
                               for b in ((3 * k, 3 * k + 1, 1), (3 * k, 3 * k + 2, 1))])),
    # ethanol: no symmetry at all, so the refinement is discrete and no search runs
    ('ethanol', dict(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])),
]


@pytest.mark.parametrize('name,case', _INVARIANCE_CASES, ids=[c[0] for c in _INVARIANCE_CASES])
def test_canonical_order_does_not_depend_on_slot_order(name, case):
    m, _ = _mol(**case)
    rng = Random(20260901)
    strings = {_canonical_string(_relabel(m, rng)) for _ in range(30)}
    assert len(strings) == 1, f'{name}: {len(strings)} distinct canonical forms'


def test_canonical_order_separates_different_structures():
    # invariance is half the contract; the other half is that the encoding still distinguishes.
    # cyclohexane against 1,5-hexadiene: same formula skeleton size, different bonds.
    ring, _ = _mol(atoms='C' * 6,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1)])
    chain, _ = _mol(atoms='C' * 6,
                    bonds=[(0, 1, 2), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2)])
    assert _canonical_string(ring) != _canonical_string(chain)


def test_canonical_order_positions_are_dense_and_zero_based():
    m, sids = _mol(atoms='C' * 8,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
                          (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
                          (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)])
    order = m.canonical_order()
    assert set(order) == set(sids)
    assert sorted(order.values()) == list(range(8))


def test_canonical_order_is_keyed_by_stable_id():
    m, sids = _mol(atoms='OCO', bonds=[(0, 1, 1), (1, 2, 1)])
    assert set(m.canonical_order()) == set(sids)
    with m.edit():
        m.delete_atom(sids[0])
    order = m.canonical_order()
    assert set(order) == {sids[1], sids[2]}
    assert sorted(order.values()) == [0, 1]


def test_canonical_order_of_empty_and_single_atom():
    m = MoleculeContainer()
    assert m.canonical_order() == {}
    m, sids = _mol(atoms='C', bonds=[])
    assert m.canonical_order() == {sids[0]: 0}


def test_canonical_order_reads_the_seed():
    # Isobutane's three methyls are interchangeable, so nothing about the structure can say which
    # of their three positions any one of them takes -- with no seed, which one lands where is
    # decided by the tie between equal certificates and moves when the slots move. Marking one
    # methyl by hand makes it distinguishable, and its position must then be a function of the
    # seeded structure alone: constant across every relabeling that carries the mark along.
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    rng = Random(20260901)
    marked, plain, orbit_positions = set(), set(), set()
    for _ in range(20):
        perm = list(m.atom_numbers)
        rng.shuffle(perm)
        out = MoleculeContainer()
        with out.edit():
            fresh = {sid: out.add_atom(m.element_of(sid)) for sid in perm}
            for a, b, order in _edges(m):
                out.add_bond(fresh[a], fresh[b], order)
        seed = {fresh[sid]: (2 if sid == sids[1] else 1) for sid in sids}
        marked.add(out.canonical_order(seed=seed)[fresh[sids[1]]])
        order = out.canonical_order()
        plain.add(order[fresh[sids[1]]])
        orbit_positions.add(frozenset(order[fresh[sid]] for sid in sids[1:]))
    assert len(marked) == 1, 'the seeded position of a distinguished methyl must be fixed'
    # Unseeded, that methyl has no position of its own: it shares an orbit with the other two, so
    # only the SET of positions the three of them occupy is a function of the structure -- which of
    # the three any one methyl takes is settled by a tie between equal certificates. As observed
    # today the unseeded position does move across relabelings, but that is not asserted: making
    # the within-orbit choice deterministic would be an improvement, not a regression.
    assert plain <= {0, 1, 2, 3}
    assert len(orbit_positions) == 1, 'the positions an orbit occupies must be fixed'


def test_canonical_order_survives_many_identical_fragments():
    # 10 cyclopropanes + 10 cyclobutanes + 10 cyclopentanes, 120 atoms of ordinary chemistry.
    #
    # This is the regression guard for the node-invariant prune. The orbit prune alone leaves three
    # candidates at every level here -- one per ring size, and no automorphism relates them -- so
    # the tree is 3**(rings) and this record blew CANON_NODE_BUDGET outright. Ranking the children
    # by an invariant and keeping only the extremal ones turns it into a path. Both build orders
    # are checked because the answer may not depend on either.
    for sizes in ([3, 4, 5] * 10, [3] * 10 + [4] * 10 + [5] * 10):
        m = _cycloalkanes(sizes)
        order = m.canonical_order()
        assert len(order) == 120, sizes
        assert sorted(order.values()) == list(range(120)), sizes
    rng = Random(20260901)
    m = _cycloalkanes([3, 4, 5] * 4)
    assert len({_canonical_string(_relabel(m, rng)) for _ in range(10)}) == 1


def _srg16(kind):
    """One of the two strongly regular graphs on parameters (16, 6, 2, 2), as a 16-atom record.

    'shrikhande' is the Cayley graph of Z4xZ4 with connection set {+-(1,0), +-(0,1), +-(1,1)};
    'rook' is the 4x4 rook's graph, vertices adjacent when they share a row or a column. They are
    not isomorphic, and no amount of refinement or local invariant can tell them apart: every
    vertex sees 6 neighbours, every adjacent pair 2 common neighbours, every non-adjacent pair 2.
    """
    idx = {(a, b): 4 * a + b for a in range(4) for b in range(4)}
    edges = set()
    for a in range(4):
        for b in range(4):
            if kind == 'shrikhande':
                for da, db in ((1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)):
                    edges.add((idx[a, b], idx[(a + da) % 4, (b + db) % 4]))
            else:
                for c in range(4):
                    if c != b:
                        edges.add((idx[a, b], idx[a, c]))
                    if c != a:
                        edges.add((idx[a, b], idx[c, b]))
    return sorted({(min(i, j), max(i, j)) for i, j in edges})


def test_canonical_order_takes_the_extremum_over_surviving_leaves():
    # The guard for the leaf certificate comparison. Both prunes above it can leave two candidates
    # standing when they are neither related by a symmetry nor separated by the node invariant --
    # then two discrete labellings reach the bottom, and only comparing their certificates picks
    # the same one every time. Delete or invert that comparison and this test reports two forms.
    #
    # Reaching that state needs a shape where individualising two inequivalent atoms yields the same
    # refinement profile at every level, which is what strongly regular graphs are. These are abstract
    # graphs and not molecules -- a six-bonded carbon is not chemistry -- but the arena accepts any
    # graph, so they are reachable input, and no chemical shape found so far reaches this node.
    #
    # ONE SRG ON ITS OWN IS THE PRIMARY CASE, and it is one component, so the per-component
    # decomposition cannot take the node away from it: individualising a vertex of an SRG(16, 6, 2, 2)
    # leaves {v}, its 6 neighbours and the 9 others, which is already equitable -- every neighbour sees
    # 2 neighbours and 3 others, every other sees 2 neighbours and 4 others -- so refinement cannot
    # split it and the search branches again. The two-component records follow because they are the
    # shapes this guard was written against, and each half must still land one form when the record is
    # canonicalised in blocks.
    cage_a = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (5, 6), (6, 7), (7, 8), (8, 9), (9, 5),
              (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]                       # pentagonal prism
    cage_b = [(0, 4), (0, 5), (0, 8), (1, 2), (1, 3), (1, 6), (2, 3), (2, 7),
              (3, 5), (4, 7), (4, 8), (5, 9), (6, 8), (6, 9), (7, 9)]
    cases = {
        'the Shrikhande graph alone': (16, [(i, j, 1) for i, j in _srg16('shrikhande')]),
        "the 4x4 rook's graph alone": (16, [(i, j, 1) for i, j in _srg16('rook')]),
        'the two SRG(16, 6, 2, 2) graphs in one record':
            (32, [(i, j, 1) for i, j in _srg16('shrikhande')]
             + [(i + 16, j + 16, 1) for i, j in _srg16('rook')]),
        'two 3-regular ten-atom cages in one record':
            (20, [(i, j, 1) for i, j in cage_a] + [(i + 10, j + 10, 1) for i, j in cage_b]),
    }
    for name, (n, bonds) in cases.items():
        m, _ = _mol(atoms='C' * n, bonds=bonds)
        rng = Random(20260901)
        strings = {_canonical_string(_relabel(m, rng)) for _ in range(20)}
        assert len(strings) == 1, f'{name}: {len(strings)} distinct canonical forms'


def test_canonical_order_refuses_to_answer_past_its_budget():
    # The guard for the no-partial-answer rule. A truncated extremal search returns some labelling
    # in place of the canonical one and nothing at the call site can tell the difference, so the
    # only safe response is to raise. The shipped budget is 1,000,000 nodes, which no record small
    # enough for a test suite can reach, hence `_node_budget` -- see its docstring.
    m, sids = _mol(atoms='C' * 8,
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
                          (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
                          (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)])
    with pytest.raises(AutomorphismBudgetExceeded) as info:
        m.canonical_order(_node_budget=2)
    assert 'refinement nodes' in str(info.value)
    assert isinstance(info.value, RuntimeError), 'callers that catch RuntimeError must still catch'
    # The failure left nothing behind: the molecule answers correctly on the next call, and the
    # answer is the one it would have given had the budgeted call never happened.
    order = m.canonical_order()
    assert sorted(order.values()) == list(range(8))
    fresh, _ = _mol(atoms='C' * 8,
                    bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 0, 1),
                           (4, 5, 1), (5, 6, 1), (6, 7, 1), (7, 4, 1),
                           (0, 4, 1), (1, 5, 1), (2, 6, 1), (3, 7, 1)])
    assert _canonical_string(m) == _canonical_string(fresh)
    # A budget large enough for the tree is not a degraded mode: same answer as the default.
    assert m.canonical_order(_node_budget=1000) == order
    # ethanol needs no search at all, so no budget can starve it
    e, _ = _mol(atoms='CCO', bonds=[(0, 1, 1), (1, 2, 1)])
    assert sorted(e.canonical_order(_node_budget=1).values()) == [0, 1, 2]


def test_canonical_order_seed_conversion_is_the_shared_one():
    m, sids = _mol(atoms='CCC', bonds=[(0, 1, 1), (1, 2, 1)])
    with pytest.raises(KeyError):
        m.canonical_order(seed={sids[0]: 1, sids[1]: 1})
    with pytest.raises(OverflowError):
        m.canonical_order(seed=dict.fromkeys(sids, -1))
