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
Salt handling in reconstruct_mapping.

An ELN record commonly stores a salt as ONE molecule container: the reactive species plus a
spectator counterion as a second connected component. reconstruct_mapping therefore considers
each container both intact and split into components -- see the comment block in reconstruct.py.
Building these cases needs `union`, because a `.` in reaction SMILES yields separate containers.
"""
from chython import ReactionContainer, smiles
from chython.algorithms.mapping.reconstruct import _mixes_views


def _fuse(rxn_smi, *, reactants=(), products=()):
    """
    Parse ``rxn_smi`` and merge the listed index groups into single containers.

    Parsing as a reaction first is essential: it gives every molecule a disjoint atom-number
    space. Independently parsed molecules all number from 1, and the resulting collision makes
    canonicalize() read the overlap as a real atom mapping and blow up in compose().
    """
    r = smiles(rxn_smi)

    def merge(mols, groups):
        if not groups:
            return list(mols)
        out = []
        for group in groups:
            m = mols[group[0]]
            for i in group[1:]:
                m = m.union(mols[i])
            out.append(m)
        return out

    out = ReactionContainer(merge(r.reactants, reactants), merge(r.products, products))
    out.canonicalize()
    return out


def test_counterion_in_reactive_container_still_couples():
    # potassium alkyltrifluoroborate logged as one record. Handed the whole container, the reactor
    # carries the spectator K+ into the generated product, which then cannot equal the salt-free
    # target -- the match is only found by also trying the borate component on its own.
    r = _fuse('Brc1ccc(cc1)C.[B-](F)(F)(F)CCC.[K+]>>c1cc(ccc1C)CCC', reactants=[(0,), (1, 2)])
    assert len(r.reactants) == 2
    assert len(r.reactants[1].connected_components) == 2, 'salt must stay one container'
    assert r.reconstruct_mapping() == ['react:suzuki']


def test_counterion_retained_on_both_sides_uses_intact_container():
    # Regression guard: when the counterion survives into the product, only the INTACT container
    # reconstructs it -- the reactive component alone generates a salt-free product that cannot
    # match. Components must therefore be additional candidates, never replacements.
    r = _fuse('OCCC.Cl>>O=CCC.Cl', reactants=[(0, 1)], products=[(0, 1)])
    assert r.reconstruct_mapping() == ['oxidize:alcohol_to_aldehyde']

    r = _fuse('CC(C)(C)OC(=O)NCCC.Cl>>NCCC.Cl', reactants=[(0, 1)], products=[(0, 1)])
    assert r.reconstruct_mapping() == ['deprotect:amine_boc']


def test_salt_break_is_not_a_reaction():
    # amine.HCl -> free amine is a workup, not a transformation: the freed component equals the
    # product, so the purification precheck rejects it.
    r = _fuse('NCCC.Cl>>NCCC', reactants=[(0, 1)])
    assert r.reconstruct_mapping() == []


def test_co_formulated_record_matches_same_as_split_record():
    # A record holding BOTH coupling partners must reconstruct exactly like the two-container
    # form: components of one source may legitimately be combined with each other.
    split = _fuse('ClC(=O)C.NCC>>CCNC(=O)C')
    fused = _fuse('ClC(=O)C.NCC>>CCNC(=O)C', reactants=[(0, 1)])
    assert len(fused.reactants) == 1
    expected = split.reconstruct_mapping()
    assert expected, 'two-container control must match'
    assert fused.reconstruct_mapping() == expected


def test_mixes_views_rejects_container_plus_own_component():
    # unit 0 = intact record A, units 1/2 = its components; unit 3 = an unrelated record B.
    a, b = object(), object()
    sources = [a, a, a, b]
    intact = [True, False, False, True]

    # intact A with either of its own components -> same atoms twice
    assert _mixes_views((0, 1), sources, intact)
    assert _mixes_views((0, 2), sources, intact)
    # two components of A -> genuine co-formulated pair
    assert not _mixes_views((1, 2), sources, intact)
    # anything from A with unrelated B
    assert not _mixes_views((0, 3), sources, intact)
    assert not _mixes_views((1, 3), sources, intact)
