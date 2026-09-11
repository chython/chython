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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""Freeze `STRUCT_VERSION == 3` arena bytes, and the answers a v3 build gave for them.

DO NOT RUN THIS AGAINST A v4 BUILD. The fixtures it writes are the evidence that a v4 reader
handles buffers it did not write; regenerating them on a v4 build would replace real v3 bytes with
v4 bytes and the compatibility suite would then prove nothing while still passing. The output is
committed for exactly that reason -- it is a frozen artifact, not a build product.

Every record is captured TWICE:

  `cold` -- serialised before any read, and
  `warm` -- serialised after `stereo_units()`, `component_labels()` and `canonical_order()`.

Those two differ on v3, for all seven records, which is the v3 defect this branch fixes: a read
builds a lazy segment, writes its table entry and bumps `total_len`, and all of that lives inside the
persistent prefix that `to_bytes()` returns. Both forms are in the wild, so a v4 reader must accept
both and give the same molecule -- in particular it must IGNORE the derived table entries a warm v3
buffer carries, which name offsets that are live-looking and meaningless.

Run from the repository root:  python chython/core/test/gen_v3_fixtures.py
"""
import base64
import textwrap

from chython.core import MoleculeContainer


def _mol(build):
    m = MoleculeContainer()
    with m.edit():
        build(m)
    return m


def ethanol_with_a_parity(m):
    # CC(O)Cl-ish: a real tetrahedral centre. C1-C2(-O3)(-Cl4) plus an H on C2.
    a = m.add_atom(6, implicit_h=3)
    b = m.add_atom(6, implicit_h=1)
    o = m.add_atom(8, implicit_h=1)
    cl = m.add_atom(17, implicit_h=0)
    m.add_bond(a, b)
    m.add_bond(b, o)
    m.add_bond(b, cl)
    m.set_parity(b, 1)


def with_coordinates(m):
    a = m.add_atom(6, implicit_h=3)
    b = m.add_atom(8, implicit_h=1)
    m.add_bond(a, b)
    m.set_xy(a, 0.0, 0.0)
    m.set_xy(b, 1.54, 0.0)


def with_stereo_groups(m):
    a = m.add_atom(6, implicit_h=3)
    b = m.add_atom(6, implicit_h=1)
    o = m.add_atom(8, implicit_h=1)
    cl = m.add_atom(17)
    m.add_bond(a, b)
    m.add_bond(b, o)
    m.add_bond(b, cl)
    m.set_parity(b, 1)
    m.set_stereo_group(b, 3, 1)      # AND1 -- racemic


def with_everything(m):
    # 2,3-dichlorobutane: two centres, coordinates, wedges, an AND group and an ABS group.
    c1 = m.add_atom(6, implicit_h=3)
    c2 = m.add_atom(6, implicit_h=1)
    c3 = m.add_atom(6, implicit_h=1)
    c4 = m.add_atom(6, implicit_h=3)
    x2 = m.add_atom(17)
    x3 = m.add_atom(17)
    m.add_bond(c1, c2)
    m.add_bond(c2, c3)
    m.add_bond(c3, c4)
    m.add_bond(c2, x2)
    m.add_bond(c3, x3)
    for i, sid in enumerate((c1, c2, c3, c4, x2, x3)):
        m.set_xy(sid, i * 1.2, (i % 2) * 0.7)
    m.set_parity(c2, 1)
    m.set_parity(c3, 2)
    m.set_wedge(c2, x2, 1)
    m.set_stereo_group(c2, 3, 1)
    m.set_stereo_group(c3, 1, 0)


def a_salt(m):
    # sodium acetate: two components, one of them a lone atom.
    c1 = m.add_atom(6, implicit_h=3)
    c2 = m.add_atom(6)
    o1 = m.add_atom(8)
    o2 = m.add_atom(8, charge=-1)
    na = m.add_atom(11, charge=1)
    m.add_bond(c1, c2)
    m.add_bond(c2, o1, 2)
    m.add_bond(c2, o2)
    del na


def a_ring(m):
    # naphthalene, Kekule, so the ring bitmap and relevant rings are non-trivial.
    ring = [m.add_atom(6, implicit_h=0) for _ in range(10)]
    bonds = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1),
             (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 0, 1), (4, 9, 1)]
    for a, b, o in bonds:
        m.add_bond(ring[a], ring[b], o)
    hs = {0: 1, 1: 1, 2: 1, 3: 1, 5: 1, 6: 1, 7: 1, 8: 1}
    for i, h in hs.items():
        m.set_hydrogens(ring[i], h)


def a_bigger_one(m):
    # a 60-atom polymethylcyclo chain, connected, with parities on every third carbon.
    prev = None
    first = None
    marks = []
    for i in range(20):
        c = m.add_atom(6, implicit_h=1)
        me = m.add_atom(6, implicit_h=3)
        m.add_bond(c, me)
        ch2 = m.add_atom(6, implicit_h=2)
        m.add_bond(c, ch2)
        if prev is not None:
            m.add_bond(prev, c)
        else:
            first = c
        prev = ch2
        marks.append(c)
    m.add_bond(prev, first)
    for c in marks:
        m.set_parity(c, 1)


CASES = [
    ('ethanol_with_a_parity', ethanol_with_a_parity),
    ('with_coordinates', with_coordinates),
    ('with_stereo_groups', with_stereo_groups),
    ('with_everything', with_everything),
    ('a_salt', a_salt),
    ('a_ring', a_ring),
    ('a_bigger_one', a_bigger_one),
]


def snapshot(m):
    """Every answer the compat test will re-check, computed by the v3 build."""
    numbers = list(m.atom_numbers)
    return {
        'atom_count': m.atom_count,
        'bond_count': m.bond_count,
        # the frozen fixtures' key; NOT an identifier -- see v3_fixtures.py's header
        'stable_ids': numbers,
        'union_feature_words': m._union_feature_words,
        'elements': [m.element_of(s) for s in numbers],
        'charges': [m.charge_of(s) for s in numbers],
        'implicit_h': [m.implicit_h_of(s) for s in numbers],
        'explicit_h': [m.explicit_h_of(s) for s in numbers],
        'hybridization': [m.hybridization_of(s) for s in numbers],
        'heteroatoms': [m.heteroatoms_of(s) for s in numbers],
        'degree': [m.degree_of(s) for s in numbers],
        'parity': [m.parity_of(s) for s in numbers],
        'stereo': [m.stereo_of(s) for s in numbers],
        'in_ring': [m.in_ring_of(s) for s in numbers],
        'ring_sizes': [m.ring_sizes_of(s) for s in numbers],
        'has_coordinates': m.has_coordinates,
        'xy': [m.xy_of(s) for s in numbers] if m.has_coordinates else None,
        'has_stereo_groups': m.has_stereo_groups,
        'stereo_groups': [m.stereo_group_of(s) for s in numbers],
        'wedges': sorted(m.wedges()),
        'bonds': sorted((min(a, b), max(a, b), m.order_of(a, b))
                        for a in numbers for b in m.neighbors_of(a)),
        'components': m.connected_components_count,
        'component_labels': m.component_labels(),
        'rings_count': m.rings_count,
        'sssr': sorted(tuple(r) for r in m.sssr),
        'atoms_order': m.atoms_order,
        'canonical_order': m.canonical_order(),
        'stereo_units': sorted(u['anchor'] for u in m.stereo_units()),
        'stereogenic_units': sorted(u['anchor'] for u in m.stereogenic_units()),
        'chiral_atoms': sorted(m.chiral_atoms()),
        'chiral_bonds': sorted(m.chiral_bonds()),
        'stereo_truncated': m.stereo_truncated,
        'validate_stereo': m.validate_stereo(),
        'features_first': m.features_of(numbers[0]),
    }


def main():
    out = []
    out.append('# Frozen v3 (`STRUCT_VERSION == 3`) arena bytes, and the answers the v3 build gave')
    out.append('# for each. GENERATED by gen_v3_fixtures.py against the last v3 build (75f6c60);')
    out.append('# never regenerate against a v4 build -- the whole point is that these bytes')
    out.append('# predate the format change.')
    out.append('')
    out.append('from base64 import b64decode')
    out.append('')
    out.append('')
    out.append('V3_FIXTURES = {}')
    out.append('')
    for name, build in CASES:
        m = _mol(build)
        cold = m.to_bytes()
        # A read builds derived segments, writes their table entries and bumps total_len -- all
        # inside the persistent prefix. So `warm` is the SAME molecule serialised differently.
        # This is the v3 defect item 3 fixes, and it is captured here on purpose: a v4 reader must
        # accept both, because both are in the wild.
        m.stereo_units()
        m.component_labels()
        m.canonical_order()
        warm = m.to_bytes()
        snap = snapshot(_mol(build))
        assert m.to_bytes()[128:] == cold[128:], name
        out.append('V3_FIXTURES[%r] = {' % name)
        out.append("    'cold': %s," % _b64(cold))
        out.append("    'warm': %s," % _b64(warm))
        out.append("    'answers': %s," % repr(snap))
        out.append('}')
        out.append('')
        print(name, len(cold), len(warm), 'differ' if cold != warm else 'IDENTICAL')
    with open('chython/core/test/v3_fixtures.py', 'w') as f:
        f.write('\n'.join(out))


def _b64(data):
    text = base64.b64encode(data).decode()
    lines = textwrap.wrap(text, 88)
    return "b64decode(\n" + '\n'.join("        %r" % chunk for chunk in lines) + '\n    )'


if __name__ == '__main__':
    main()
