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
"""Every `z` translation in `Z3_MAP`, checked against a representative molecule per pattern.

chython 2 saturates `z` and caps it at 3; the core reports what it found -- 3 is sp and nothing else,
5 is two cumulated doubles with no triple, 6 is any other combination.  Copying a `z3` across
narrows the pattern, and a narrowed prefilter does not raise: the rule silently stops firing.
"""
import pytest

from .gen_standardize_rules import Z3_MAP
from .._smarts import compile_smarts
from chython.core import read_smiles


# V2 pattern -> a molecule the V2 pattern matches.  `+` marks the seven written by hand
# because nothing in the standardization corpus reaches that pattern.
REPRESENTATIVES = {
    # Group A -- genuine sp, `z3` stays `z3`: every one of these writes a triple bond on the `z3` atom
    '[N;D2;z3;+](#[N;D1])[C,N,O;z1;-]': '[CH2-][N+]#N',
    '[N;D2;z3;x2]([N;D2;z1])#[N;D1]': 'CN[N]#N',
    '[N;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]': 'C#[N]O',
    '[N;D2;z3;x1]([N;D2;z1])#[C;D1,D2]': 'C#[N]NC',
    '[N;D2;z3;x1;+]([N,O,S;D1])#[C;D1;-]': '[C-]#[N+]O',
    '[N;D2;z3;x1;+]([N;D2;z1])#[C;D1;-]': '[C-]#[N+]NC',
    '[N;D2;z3]([A])#[C;D1]': 'C#[N]O',
    '[N;D3;z3;x1](#[N;D1])(C)C': 'CN(#N)C',
    '[N;D1;x0;z3]#[C;D2;z3;x2][O;D1]': 'C(#N)O',
    '[N;D1;x0;z3]#[C;D2;z3;x2][O;D1;-]': 'C(#N)[O-]',
    '[C;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]': 'C(#C)O',
    '[C;D2;z3;x1]([N;D2;z1])#[C;D1,D2]': 'C(#C)NC',

    # Group C -- two cumulated doubles, no triple.  `z3` -> `z5`
    '[N;D3;z3;x2](=[O;D1])([O;D1])=C': 'C=N(=O)O',
    '[N;D3;z3](=[O;D1])(=[C,N,O])-[A]': 'C=N(=O)O',
    '[N;D3;z3](=[N;D3;z2;+])(=[O;D1])[A]': 'CN(=O)=[N+](C)C',                            # +
    '[N;D3;z3](=[N;D3;z2;+])(=[N;D1,D2;z2])[A]': 'CN(=N)=[N+](C)C',                      # +
    '[N;D3;z3](=[N;D1,D2;z2])(=[C,N])[A]': 'CN(=N)=N',
    '[N;D3;z3](=[O;D1])(=[O;D1])[A;-]': 'N([O-])(=O)=O',
    '[N;D2;z3;x2;-](=[O;D1])=[O;D1]': '[N-](=O)=O',
    '[N;D2;z3;x2](=[N;D2;z2])=[N;D1]': 'CN=[N]=N',
    '[N;D2;z3;x1;+](=[N;D1])=[C;D1,D2;z2;-]': '[CH-]=[N+]=N',
    '[P;D4;z3;-](=[O;D1])(=[O;D1])([A])[A]': 'C[P-](=O)(=O)C',                           # +
    '[S;D1;-][S;D4;z3](=[O;D1])(=[A])[A]': '[S-]S(=O)(=C)C',                             # +
    '[S;D1][S;D4;z3](=[O;D1])(=[A])[A]': 'SS(=O)(=C)C',                                  # +
    '[S;D3;z3;-](=[O;D1])(=[O;D1])[A]': 'C[S-](=O)=O',
    '[S;D3;z3;x3;-]([S;D1;-])(=[O;D1])=[O;D1]': '[S-][S-](=O)=O',                        # +
    '[S;D4;z3:1]([O;D1:2])(=[N;D1,D2;z2:3])(=[A])[A]': 'CN=S(=N)(O)O',

    # Group D -- anything else the core did not fold into 5.  `z3` -> `z6`
    '[N;D2;z3](#[N;D1])=[C,N,O]': 'C=N#N',
    '[N;D2;z3;x2](#[N;D2;+][A])=[N;D1;-]': 'C[N+]#N=[N-]',
    '[N;D2;z3;x2](=[N;D2;z2])#[N;D1;-]': 'CN=N#[N-]',
    '[N;D2;z3;x2](=[N;D1;-])#[N;D1]': '[N-]=N#N',
    '[N;D2;z3;x1](=[N;D1])#[C;D1,D2]': 'C#N=N',
    '[N;D2;z3;x1](=[N,O;z2])#[C;D1,D2]': 'C#N=N',
    '[S;D4;z3;-](=[O;D1])(=[O;D1])(=[O;D1])[A]': 'C[S-](=O)(=O)=O',                      # +
}


def _sites(pattern: str, source: str):
    """Which atom sets the pattern covers in the core.  A set, because automorphisms repeat one."""
    query, _, _ = compile_smarts(pattern)
    molecule = read_smiles(source)
    return {tuple(sorted(mapping.values())) for mapping in query.get_mapping(molecule)}


def test_every_translated_primitive_has_a_representative():
    """`Z3_MAP` and `REPRESENTATIVES` cover each other exactly, in both directions."""
    assert set(REPRESENTATIVES) == set(Z3_MAP), (
        'Z3_MAP and REPRESENTATIVES disagree.\n  no representative: '
        f'{sorted(set(Z3_MAP) - set(REPRESENTATIVES))}\n  no ruling: '
        f'{sorted(set(REPRESENTATIVES) - set(Z3_MAP))}')


@pytest.mark.parametrize('pattern', sorted(p for p, (g, _) in Z3_MAP.items() if g != 'A'))
def test_the_untranslated_pattern_would_have_dropped_the_match(pattern):
    """An untranslated `z3` matches nothing, so the translation was necessary."""
    source = REPRESENTATIVES[pattern]
    assert not _sites(pattern, source), (
        f'{pattern!r} still matches {source!r} in the core with its `z3` untouched.  Either the '
        "core's `z3` has widened -- in which case Z3_MAP needs re-deriving, not this test relaxing "
        '-- or this representative no longer isolates the primitive')


def test_group_a_keeps_z3_because_a_triple_bond_is_written_beside_it():
    """Group A keeps `z3` only because every one of its patterns writes an explicit `#` bond.

    The core's `z3` is the narrowest of the three values, so leaving one alone is the riskiest
    outcome; the explicit triple already pinned the atom to sp in chython 2 too.
    """
    for pattern, (group, _) in Z3_MAP.items():
        if group == 'A':
            assert '#' in pattern, (
                f'{pattern!r} keeps `z3` without an explicit triple bond to justify it.  The core '
                "reads `z3` as sp and nothing else, so this needs measuring, not inheriting")
