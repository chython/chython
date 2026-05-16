# -*- coding: utf-8 -*-
#
#  Copyright 2026 Tagir Akhmetshin <tagirshin@gmail.com>
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
"""Regression tests for the SMARTS writer's handling of CXSMILES extension
blocks (`|...|`).

Per the CXSMILES spec, a CXSMILES extension block, if present, must appear
**at the very end** of the SMILES/SMARTS string, after a single space — never
between disconnected fragments, never twice in one string, and never inside a
reaction string (i.e. not before a `>`).
"""
import re

from chython import smarts, ReactionContainer
from chython.containers.query import QueryContainer
from chython.periodictable import QueryElement


_CX_RE = re.compile(r"\|[^|]*\|")


def _assert_cxsmiles_well_formed(out: str):
    """Generic invariant: any `|...|` block in ``out`` must be the trailing
    token, preceded by a single space, and there must be at most one such
    block in the whole string."""
    blocks = list(_CX_RE.finditer(out))
    if not blocks:
        return
    assert len(blocks) == 1, (
        f"CXSMILES blocks must be consolidated; found multiple in {out!r}"
    )
    block = blocks[0]
    assert block.end() == len(out), (
        f"CXSMILES block must be at the end of the string, got {out!r}"
    )
    assert out[block.start() - 1] == " ", (
        f"CXSMILES block must be preceded by a single space, got {out!r}"
    )
    # the part before the block must contain no '>' followed by anything
    # (i.e. the CX must trail the whole reaction, not a single side)
    if ">" in out:
        # everything between the LAST '>' and the block must be the products side
        last_gt = out.rfind(">")
        assert last_gt < block.start(), (
            f"CXSMILES block sits before a reaction separator in {out!r}"
        )


def _make_radical_query(atom_n: int, mapping: int) -> QueryContainer:
    """Build a 2-atom query whose second atom carries a radical mark."""
    q = QueryContainer("")
    q.add_atom(QueryElement.from_atomic_number(atom_n)(), mapping)
    e = QueryElement.from_atomic_number(atom_n)()
    e._is_radical = True
    q.add_atom(e, mapping + 1)
    q.add_bond(mapping, mapping + 1, 1)
    q._smarts = ""  # force regeneration (str() prefers cached _smarts)
    return q


def test_cxsmiles_block_only_trails_on_single_fragment_smarts():
    # Single-fragment query with a radical: the |^1:N| extension must trail
    # the SMARTS body, after exactly one space.
    q = _make_radical_query(7, 1)
    out = str(q)
    _assert_cxsmiles_well_formed(out)
    # Sanity: there IS a CX block (radical present).
    assert "|" in out


def test_cxsmiles_block_consolidated_across_fragments_in_one_side():
    # Two disconnected fragments where the *second* fragment has a radical.
    # The CXSMILES block must not appear between the fragments.
    qA = QueryContainer("")
    qA.add_atom(QueryElement.from_atomic_number(6)(), 1)
    qA.add_atom(QueryElement.from_atomic_number(6)(), 2)
    qA.add_bond(1, 2, 1)
    qA._smarts = ""

    qB = _make_radical_query(7, 3)

    rxn = ReactionContainer([qA, qB], [qA])
    out = str(rxn)
    _assert_cxsmiles_well_formed(out)


def test_cxsmiles_block_consolidated_across_reaction_sides():
    # Reaction whose product side has a radical. The CXSMILES block must
    # trail the *whole* reaction, not appear mid-side.
    qA = QueryContainer("")
    qA.add_atom(QueryElement.from_atomic_number(6)(), 1)
    qA.add_atom(QueryElement.from_atomic_number(6)(), 2)
    qA.add_bond(1, 2, 1)
    qA._smarts = ""

    # product side: same C-C plus a radical N-O fragment.
    qB = _make_radical_query(7, 3)
    qC = QueryContainer("")
    qC.add_atom(QueryElement.from_atomic_number(6)(), 5)
    qC.add_atom(QueryElement.from_atomic_number(6)(), 6)
    qC.add_bond(5, 6, 1)
    qC._smarts = ""

    rxn = ReactionContainer([qA], [qC, qB])
    out = str(rxn)
    _assert_cxsmiles_well_formed(out)


def test_cxsmiles_roundtrip_buggy_input_normalized():
    # The buggy TSV row from the priority-rules report had the CXSMILES
    # block placed mid-fragment (between disconnected pieces of one
    # side). Whatever we read in, the writer must normalize the output
    # to put the block at a single trailing position.
    body = (
        "[C;D3:3](=[C;D2:4])(-[C;D3:5])-[Cl;D1:6] |^1:2|"
        ".[C;D2:1]-[N;D1:2]"
    )
    rxn = smarts(f"{body}>>{body}")
    out = str(rxn)
    _assert_cxsmiles_well_formed(out)


def test_cxsmiles_absent_when_no_radicals():
    # A query without any radical/stereo-extension content must produce
    # no CXSMILES extension at all.
    q = QueryContainer("")
    q.add_atom(QueryElement.from_atomic_number(6)(), 1)
    q.add_atom(QueryElement.from_atomic_number(6)(), 2)
    q.add_bond(1, 2, 1)
    q._smarts = ""
    out = str(q)
    assert "|" not in out


def test_cxsmiles_indices_reference_correct_atoms():
    # The radical indices in the consolidated CX block must point at
    # *radical* atoms in the emitted SMARTS — round-tripping through
    # smarts() must preserve the set of radical atoms.
    qA = QueryContainer("")
    qA.add_atom(QueryElement.from_atomic_number(6)(), 1)
    qA.add_atom(QueryElement.from_atomic_number(6)(), 2)
    qA.add_bond(1, 2, 1)
    qA._smarts = ""

    qB = _make_radical_query(7, 3)
    rxn = ReactionContainer([qA, qB], [qA])
    out = str(rxn)
    rxn2 = smarts(out)
    radicals_in = sum(
        1 for m in rxn.reactants for _, a in m.atoms() if getattr(a, "is_radical", False)
    )
    radicals_out = sum(
        1 for m in rxn2.reactants for _, a in m.atoms() if getattr(a, "is_radical", False)
    )
    assert radicals_in == radicals_out, (
        f"radical count must roundtrip, in={radicals_in}, out={radicals_out}"
    )
