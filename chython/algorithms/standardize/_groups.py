# -*- coding: utf-8 -*-
#
#  Copyright 2021-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy


def _rules():
    """
    rules without overlapping. these rules can match once to same set of atoms.
    """
    from ... import smarts
    from ...containers import QueryContainer

    rules = []

    #
    #     A         A
    #     |         |
    #     P    >>  [P+]
    #   / | \     / | \
    #  A  A  A   A  A  A
    #
    q = smarts('[P;D4;x0;z1]')
    atom_fix = {1: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # A   H   A     A     H   A
    #  \ / \ /       \  .. \ /
    #   B   B    >>   B     B
    #  / \ / \       / \  .. \
    # A   H   A     A   H     A
    #
    q = smarts('[B;z1:1]1[H;D2:3][B;z1:2][H;D2:4]1')
    atom_fix = {}
    bonds_fix = ((1, 3, 8), (2, 4, 8))
    rules.append((q, atom_fix, bonds_fix, False))

    # OGB DS
    #
    # O*  A  [O-]  A
    #  \ /     \  /
    #   N  >>  [N+]
    #   |       ||
    #  C,N*     C,N
    #
    q = smarts('[O;D1;z1][N;D3;z1][C,N;z1] |^1:0,2|')
    atom_fix = {1: (-1, False), 2: (1, None), 3: (0, False)}
    bonds_fix = ((2, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #       |            |
    # O* -- S -- O* >> O=S=O
    #       |            |
    #
    q = smarts('[O,S;D1;z1][S;D4;z1][O;D1;z1] |^1:0,2|')
    atom_fix = {1: (0, False), 3: (0, False)}
    bonds_fix = ((1, 2, 2), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #      A         A
    #     //        //
    # B - N >> B .. N
    #     \          \
    #      A          A
    #
    q = smarts('[B]-[N;D3;z2]')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     A             A
    #     |             |
    # B = N - A >> B .. N - A
    #     |             |
    #     A             A
    #
    q = smarts('[B;z1,z2]-,=[N;D4;z1,z2]')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # R2S - B >> R2S .. B
    #
    q = smarts('[B;z1]-[O,S;D3;z1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  [B-] = [N+] >>  B - N
    #
    q = smarts('[B;z2;-]=[N;D1,D2,D3;z2;+]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #        [A-]                 A
    #         |                   |
    # [A-] - [B+3] - [A-] >> A - [B-] - A
    #         |                   |
    #        [A-]                 A
    #
    q = smarts('[B;D4;z1;+3]([A;-])([A;-])([A;-])[A;-]')
    atom_fix = {1: (-4, None), 2: (1, None), 3: (1, None), 4: (1, None), 5: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #        A             A
    #        |             |
    # [A-] - B - A >> A - [B-] - A
    #        |             |
    #        A             A
    #
    q = smarts('[B;D4;z1]-[A;-]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #      A             A
    #      |             |
    # A -  B - A >> A - [B-] - A
    #      |             |
    #      A             A
    #
    q = smarts('[B;D4;z1]')
    atom_fix = {1: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     A         A
    #     |         |
    #     N    >>  [N+]
    #   / | \     / | \
    #  A  A  A   A  A  A
    #
    q = smarts('[N;D4;z1]')
    atom_fix = {1: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #  aci-nitro
    #      O          O
    #     //         //
    # C = N  >> C - [N+]
    #      \         \
    #       OH       [O-]
    #
    q = smarts('[N;D3;z3;x2](=[O;D1])([O;D1])=C')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 4, 1),)
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #            O               [O-]
    #           //               /
    # [C,N,O] = N  >> [C,N,O] = [N+]
    #           \                \
    #            A                A
    #
    q = smarts('[N;D3;z3](=[O;D1])(=[C,N,O])-[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    # fix aci-nitro
    #       [O-]          [O-]
    #       /             /
    # C = [N+ ] >>  C - [N+]
    #      \             \\
    #       OH            O
    #
    q = smarts('[N;D3;z2;x2;+]([O;D1;-])([O;D1])=C')
    atom_fix = {}
    bonds_fix = ((1, 3, 2), (1, 4, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    # fix CN(=O)=N(=O)C
    #
    # [N+] = N = O >> [N+] = [N+] - O-
    #        |                |
    #        A                A
    #
    q = smarts('[N;D3;z3](=[N;D3;z2;+])(=[O;D1])[A]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 3, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    # Fix CN(=O)=N(=N)C
    # [N+] = N = N - ? >> [N+] = [N+] - [N-] - ?
    #        |                    |
    #        A                    A
    #
    q = smarts('[N;D3;z3](=[N;D3;z2;+])(=[N;D1,D2;z2])[A]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 3, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    # For N-case is not unique!
    #          N             [N-]
    #         //             /
    # [C,N] = N  >> [C,N] = [N+]
    #         \              \
    #          A              A
    #
    q = smarts('[N;D3;z3](=[N;D1,D2;z2])(=[C,N])[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [N-] - [N+] = O >> N = [N+] - [O-]
    #         |               |
    #         A               A
    #
    q = smarts('[N;D3;z2;+](=[O;D1])[N;D1,D2;z1;-]')
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #         O              [O-]
    #        //              /
    # [A-] - N   >> [A-] - [N+]
    #        \\             \\
    #         O              O
    #
    q = smarts('[N;D3;z3](=[O;D1])(=[O;D1])[A;-]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # O : N : O      O = [N+] - [O-]
    #     |      >>       |
    #     A               A
    #
    q = smarts('[N;D3;a](:[O;D1])(:[O;D1])[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    # Nitrite
    #
    #   O        [O-]
    #  //        /
    # [N-]  >>  N
    #  \\       \\
    #   O        O
    #
    q = smarts('[N;D2;z3;x2;-](=[O;D1])=[O;D1]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [O,C,N] = N # N >> [O,C,N] = [N+] = [N-]
    #
    q = smarts('[N;D2;z3](#[N;D1])=[C,N,O]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [C,N,O-] - [N+] # N  >> [C,N,O] = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;+](#[N;D1])[C,N,O;z1;-]')
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  A - [N+] # N = [N-]  >> A - N = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;x2](#[N;D2;+][A])=[N;D1;-]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # A - N = N = N >> A - N = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;x2](=[N;D2;z2])=[N;D1]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # A - NH - N # N >> A - N = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;x2]([N;D2;z1])#[N;D1]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [N-] # N = N - A >> [N-] == [N+] == N - A
    #
    q = smarts('[N;D2;z3;x2](=[N;D2;z2])#[N;D1;-]')
    atom_fix = {1: (1, None)}
    bonds_fix = ((1, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [N-] == N # N >> [N-] == [N+] == [N-]
    #
    q = smarts('[N;D2;z3;x2](=[N;D1;-])#[N;D1]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # A - C # N = NH >> A - [CH] = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;x1](=[N;D1])#[C;D1,D2]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    # note: order dependent
    # A - C # N = [O,N] >> A - C # [N+] - [O,N-]
    #
    q = smarts('[N;D2;z3;x1](=[N,O;z2])#[C;D1,D2]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [NH2,OH,SH] - N # C - A >> [NH,O,S-] - [N+] # C - A
    #
    q = smarts('[N;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # A - [NH] - N # C >> A - [N-] - [N+] # C
    #
    q = smarts('[N;D2;z3;x1]([N;D2;z1])#[C;D1,D2]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [NH2,OH,SH] - [N+] # [C-] >> [NH,O,S-] - [N+] # [CH]
    #
    q = smarts('[N;D2;z3;x1;+]([N,O,S;D1])#[C;D1;-]')
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # A - [NH] - [N+] # [C-] >> A - [N-] - [N+] # [CH]
    #
    q = smarts('[N;D2;z3;x1;+]([N;D2;z1])#[C;D1;-]')
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # A - N # C >> A - [N+] # [C-]
    #
    q = smarts('[N;D2;z3]([A])#[C;D1]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    # fix old diazo rule
    #
    # A - [C-] = [N+] = [NH] >> A - [CH] = [N+] = [N-]
    #
    q = smarts('[N;D2;z3;x1;+](=[N;D1])=[C;D1,D2;z2;-]')
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #    |           |
    #  - N -   >> - [N+] -
    #    \\          |
    #    [O,N]      [O,N-]
    #
    q = smarts('[N;D4;z2]=[O,N;z2]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    # N-oxide radical
    #
    #    |         |
    #  - N*  >>  - N
    #    \\        |
    #     O        O*
    #
    q = smarts('[N;D3;z2]=[O;D1] |^1:0|')
    atom_fix = {1: (0, False), 2: (0, True)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    # false N-oxide radical
    #
    #    |        |
    #  = N  >> = [N+]
    #     \        \
    #      O*       [O-]
    #
    q = smarts('[O;D1][N;D3;z2] |^1:0|')
    atom_fix = {1: (-1, False), 2: (1, False)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # C           C
    #  \           \
    #   N # N >>   [N+] = [N-]
    #  /           /
    # C           C
    #
    q = smarts('[N;D3;z3;x1](#[N;D1])(C)C')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  C - N = [N+]  >>  C - [N+] # N
    #
    q = smarts('[N;D1;z2;x1;+]=[N;D2;x1;z2]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 3),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  [N+] - [C-] = O  >>  N = C = O
    #
    q = smarts('[C;D2;z2;x2;-]([N;D1,D2;z1;+])=[O;D1]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  N # C - OH  >>  HN = C = O
    #
    q = smarts('[N;D1;x0;z3]#[C;D2;z3;x2][O;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #  N # C - [O-]  >>  [N-] = C = O
    #
    q = smarts('[N;D1;x0;z3]#[C;D2;z3;x2][O;D1;-]')
    atom_fix = {1: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # - [N+] - [O-]  >>  - N = O
    #
    q = smarts('[O;D1;z1;x1;-][N;D2;z1;+]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # CH-N=O >> C=N-OH
    #
    q = smarts('[O;D1;z2;x1]=[N;D2;x1;z2][C;D1,D2,D3;z1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #       [O+]R           [O]R
    #       //              /
    # R2N - C  >> R2[N+] = C
    #
    q = smarts('[O;D2;z2;+]=[C;z2][N;D1,D2,D3;z1]')
    atom_fix = {1: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # fix pyrylium dearomatization. restore O+ in aromatic ring from N+ form.
    # e.g. O1C(=[NH+]...)... -> [O+]1=C(N...)...
    #
    q = smarts('[O;D2;r6]1[C;z2](=[N;+])[A;z2]-,=[A;z2]-,=[A;z2]-,=[A;z2]1')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #      [O,S]H    [O,S]
    #      /         //
    # N = C  >> NH - C
    #
    q = smarts('[N;z2]=[C;D2,D3;z2]-[O,S;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #       OH          O
    #      /           //
    # C = C    >> C - C
    #      \           \
    #      [O,N]       [O,N]
    #
    q = smarts('[O;D1;x0;z1]-[C;D3;z2;x2](-[O,N])=C')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 4, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix 1,3,5-triketone (cyanuric acid type)
    # e.g. O=C1NC(=O)NC(=O)N1 → Oc1nc(O)nc(O)n1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2]C(=[O,S,N;D1])[C,N;h1,h2]C(=[O,S,N;D1])[C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (4, 5, 1), (4, 6, 2), (7, 8, 1), (7, 9, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-1,3-dione (uracil-type)
    # e.g. O=C1NC(=O)C=CN1 → Oc1nc(O)ccn1 (both keto→enol, restores aromaticity)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2]C(=[O,S,N;D1])[A;z2]-,=[A;z2][C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 8, 2), (4, 5, 1), (3, 4, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-1,2-dione
    # e.g. O=C1C(=O)NNC=C1 → Oc1c(O)nncc1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1C(=[O,S,N;D1])[C,N;h1,h2][A;z2]-,=[A;z2][C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 8, 2), (3, 4, 1), (3, 5, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-1,4-dione
    # e.g. O=C1NNC(=O)C=C1 → Oc1nnc(O)cc1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][C,N;h1,h2]C(=[O,S,N;D1])[A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (5, 6, 1), (4, 5, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-2-one. note: only after amide rule
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2;z1][A;z2]-,=[A;z2][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-4-one
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[A;z2]=[A;z2][N;D2;z1][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-2,5-dione
    # e.g. O=C1CCC(=O)N1 → OC1=CC=C(O)N1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[N,O,S;z1;D2,D3]C(=[O,S,N;D1])[C,N;h1,h2][C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 7, 2), (4, 5, 1), (4, 6, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-2,4-dione
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[N,O,S;z1;D2,D3][C,N;h1,h2]C(=[O,S,N;D1])[C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 7, 2), (5, 6, 1), (4, 5, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-2-one
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N,O,S;z1;D2,D3][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-1-one
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][A;z2]=,-[A;z2]-[N,O,S;z1;D2,D3]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-2-imine. move double bond to the ring
    #
    q = smarts('[N;z2;!R]=[C;D3;r6]1[C,N;h1,h2][A]-,=[A]-,=[A]-,=[A]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-2-imin. move double bond to the ring
    #
    q = smarts('[N;z2;!R]=[C;D3;r5]1[C,N;h1,h2][N,O,S;z1;D2,D3][A]-,=[A]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyrrole-2-imin. move double bond to the ring
    #
    q = smarts('[N;z2;!R]=[C;D3;r5]1[C,N;h1,h2;z1][A]-,=[A][N,O,S;z1;D2,D3]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 pyridinone. bridge N(D3,z1) adjacent to z1 atom.
    # e.g. OC1=NC2=CC=CN2C=C1 (pyrrolo[1,2-a]pyrimidinone)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][C;D3;z2;r5]=2[N;D3]([A;z2]-,=[A;z2][A;z2]=2)[A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 pyridinone. bridge at ring closure, N(D3,z1) last.
    # e.g. OC1=CC=CC2=CC=CN12 (indolizinone)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][A;z2]-,=[A;z2][C;D3;z2;r5]=2[N;D3;z1]1[A;z2]-,=[A;z2][A;z2]=2')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 pyridinone. bridge at ring closure, N(D3,z1) before C(D3).
    # e.g. OC1=NC=CN2C=CC=C12 (imidazo[1,2-a]pyridinone)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][A;z2]-,=[A;z2][N;D3;z1;r5]2[C;D3;z2]1=[A;z2][A;z2]-,=[A;z2]2')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 pyridinone. C=O adjacent to bridge N(D3,z1).
    # e.g. O=C1CN2C=CC=C2C=C1 -> OC1=CN2C=CC=C2C=C1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][N;D3;z1;r5]2[C;D3;z2](=[A;z2][A;z2]-,=[A;z2]2)[A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 diketone. two C=O in 6-ring, bridge N(D3,z1) between them.
    # e.g. O=C1CN2C=CC=C2C(=O)C1 -> Oc1cn2cccc2c(O)c1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][N;D3;z1;r5]2[C;D3;z2](=[A;z2][A;z2]-,=[A;z2]2)C(=[O,S,N;D1])[C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (9, 10, 1), (9, 11, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 diketone. C=O directly on bridge N, second C=O further.
    # e.g. O=C1N2C=CC=C2CC(=O)C1 -> Oc1n2cccc2cc(O)c1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[N;D3;z1;r5]2[C;D3;z2](=[A;z2][A;z2]-,=[A;z2]2)[C,N;h1,h2]C(=[O,S,N;D1])[C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 11, 2), (9, 10, 1), (8, 9, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 diketone. both C=O adjacent to bridge N.
    # e.g. O=C1N2C=CC=C2C(=O)CC1 -> Oc1n2cccc2c(O)cc1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[N;D3;z1;r5]2[C;D3;z2](=[A;z2][A;z2]-,=[A;z2]2)C(=[O,S,N;D1])[C,N;h1,h2][C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 11, 2), (8, 9, 1), (8, 10, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r6 diketone. two C=O in 6-ring with bridge N between them (1,4-dione with N in path).
    # e.g. OC1=NN2C=CC=C2N=C1O (degraded to O=C1NNC2=CC=CC2NC1=O)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][N;D3;z1;r5]2[C;D3;z2](=[A;z2][A;z2]-,=[A;z2]2)[C,N;h1,h2]C1=[O,S,N;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (10, 11, 1), (9, 10, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r5 azolone. C=O with C;z1 and S/O;z1 in same ring, bridge to second r5.
    # e.g. O=C1CSC2=CC=CN12 -> OC1=CSC2=CC=CN12 (indolizine-thiazolone)
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N,O,S;z1;D2,D3][C;D3;z2]=2[N;D3;z1]1[A;z2;r5]-,=[A;z2][A;z2]=2')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r5 azolone. C=O with heteroatom at ring closure, bridge N between rings.
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N;D3;z1]2[C;D3;z2](=[A;z2;r5][A;z2]-,=[A;z2]2)[N,O,S;z1;D2,D3]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r5 azolone. C=O with N bridge then C bridge (swapped).
    # e.g. O=C1CNN2C=CC=C12 -> OC1=CNN2C=CC=C12
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N,O,S;z1;D2,D3][N;D3;z1]2[A;z2;r5]-,=[A;z2][A;z2]=[C;D3]12')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix fused r5+r5 azolone. C=O with C bridge then N bridge.
    # e.g. O=C1CC2=CC=CN2N1 -> Oc1cn2cccc2[nH]1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][C;D3;z2]2=[A;z2][A;z2]-,=[A;z2][N;D3;z1]2[C,N;h1,h2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix charged azolone (imidazolium/oxazolium-type with C=O)
    # e.g. CN1N=[N+](C)CC1=O -> Cn1n[n+](C)cc1O, C[N+]1=NOC(=O)C1 -> C[n+]1noc(O)c1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N;z2;+]-,=[A;z2]-[N,O,S;z1;D2,D3]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix charged azolone (alternative isomer, N+ at ring closure)
    # e.g. CN1CC(=O)[N+](C)=N1 -> Cn1cc(O)[n+](C)n1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r5]1[C,N;h1,h2][N,O,S;z1;D2,D3][A;z2]-,=[N;z2;+]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix charged pyridinone (N+ in 6-ring with C=O)
    # e.g. C[N+]1=CC=CC(=O)N1 -> C[n+]1cc=cc(O)n1
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[C,N;h1,h2][N;z2;+]-,=[A;z2]-,=[A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    # acyclic keto-enol
    #         OH                  O
    #        /                   //
    # C,H - C = C - C,H >> C,H - C - C - C,H
    #
    q = smarts('[O;D1;z1;x0][C;D2,D3;z2;x1;!R]=[C;z2;x0]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #       A                   A
    #       |                   |
    #  A - [P+] - [O-]  >>  A - P = O
    #       |                   |
    #       A                   A
    #
    q = smarts('[P;D4;z1;+][O;D1;-]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #       A                   A
    #       |                   |
    #  A - [P-] - [C+]  >>  A - P = C
    #       |                   |
    #       A                   A
    #
    q = smarts('[P;D4;z1;-][C;D1,D2,D3;z1;+]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #   F   F         F   F
    #    \ /           \ /
    # F - P - F >> F - [P-] - F
    #    / \           / \
    #   F   F         F   F
    #
    q = smarts('[P;D6;z1]([F;D1])([F;D1])([F;D1])([F;D1])([F;D1])[F;D1]')
    atom_fix = {1: (-1, None)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     O               O
    #     \\              \\
    # A - [P-] - A  >> A - P - A
    #     //              /
    #     O              [O-]
    #
    q = smarts('[P;D4;z3;-](=[O;D1])(=[O;D1])([A])[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     O           [O-]
    #     \\            \
    # A -  P - A  >> A - P - A
    #     /             //
    #   [S-]            S
    #
    q = smarts('[S;D1;z1;x1;-][P;D4;z2]=[O;D1]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     O           [OH]
    #     \\            \
    # A -  P - A  >> A - P - A
    #     /             //
    #   [SH]            S
    #
    q = smarts('[S;D1;z1;x1][P;D4;z2]=[O;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #     A                A
    #     \\               \\
    # A -  S - [S-]  >> A - S - [O-]
    #     //               //
    #     O                S
    #
    q = smarts('[S;D1;-][S;D4;z3](=[O;D1])(=[A])[A]')
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     A                A
    #     \\               \
    # A -  S - [SH]  >> A - S - [OH]
    #     //               //
    #     O                S
    #
    q = smarts('[S;D1][S;D4;z3](=[O;D1])(=[A])[A]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #           A                   A
    #          /                   /
    # [O-] - [S,Si,Se+]  >>  O = [S,Si,Se]
    #          \                   \
    #           A                   A
    #
    q = smarts('[S,Se,Si;D3;z1;+][O;D1;-]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #       A?                  A?
    #       |                   |
    #  A = [S+] - [O-]  >>  A = S = O
    #       |                   |
    #       A?                  A?
    #
    q = smarts('[S;D2,D4;z2;+][O;D1;-]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #       [O-]          O
    #       /            //
    # A = [S+2]  >>  A = S
    #       \            \\
    #       [O-]          O
    #
    q = smarts('[S;D3;z2;+2]([O;D1;-])[O;D1;-]')
    atom_fix = {1: (-2, None), 2: (1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #          A                  A
    #          |                  |
    # [O-] - [S+2] - [O-] >>  O = S = O
    #          |                  |
    #          A                  A
    #
    q = smarts('[S;D4;z1;+2]([O;D1;-])[O;D1;-]')
    atom_fix = {1: (-2, None), 2: (1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  A - [S-] - [C+]  >>  A - S = C
    #       |                   |
    #       A                   A
    #
    q = smarts('[S;D3;z1;-]([C;D1,D2,D3;z1;+])')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  O              O
    #  \\             \\
    #  [S-] - A  >>    S - A
    #  //             /
    #  O            [O-]
    #
    q = smarts('[S;D3;z3;-](=[O;D1])(=[O;D1])[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #     O             [O-]
    #     \\             \
    # A - [S-] = O >> A - S = O
    #     //             //
    #     O              O
    #
    q = smarts('[S;D4;z3;-](=[O;D1])(=[O;D1])(=[O;D1])[A]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  O              [O-]
    #  \\              \
    #  [S-] - [S-] >>   S = S
    #  //              /
    #  O              [O-]
    #
    q = smarts('[S;D3;z3;x3;-]([S;D1;-])(=[O;D1])=[O;D1]')
    atom_fix = {1: (1, None), 2: (1, None), 3: (-1, None), 4: (-1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 1), (1, 4, 1))
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #    A            A
    #    \            \
    # A - S = O >> A - S = O
    #    /            //
    #  [OH]           O
    #
    q = smarts('[S;D4;z2](=[O;D1])[O;D1]')
    atom_fix = {}
    bonds_fix = ((1, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #    A            A
    #    \            \
    # A - S = O >> A - S = O
    #    /            //
    #  [NH]           N
    #
    q = smarts('[S;D4;z2](=[O;D1])[N;D1,D2;z1]')
    atom_fix = {}
    bonds_fix = ((1, 3, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #      A{1,3}               A{1,3}
    #      |                    |
    #  N = S - [OH]  >>  [NH] - S = O
    #
    q = smarts('[S;D3,D5;z2](=[N;D1,D2;z2])[O;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    #     [OH]                O
    #      |                 //
    #  N = S = A  >>  [NH] - S = A   (duplicated for overlapping groups)
    #      |                 |
    #      A                 A
    #
    q = smarts('[S;D4;z3:1]([O;D1:2])(=[N;D1,D2;z2:3])(=[A])[A]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 1))
    rules.append((q, atom_fix, bonds_fix, True))
    rules.append((q, atom_fix, bonds_fix, True))  # second shot for overlapping groups

    #
    # C # C - [O,NH,S]H  >> C=C=[O,NH,S]
    #
    q = smarts('[C;D2;z3;x1]([N,O,S;D1])#[C;D1,D2]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # C # C - [NH]R  >> C=C=NR
    #
    q = smarts('[C;D2;z3;x1]([N;D2;z1])#[C;D1,D2]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    # Carbon Monoxide
    #
    # [CX1] = O  >> [С-] # [O+]
    #
    q = smarts('[C;D1;x1;z2]=[O;D1] |^1:0|')
    atom_fix = {1: (-1, False), 2: (1, None)}
    bonds_fix = ((1, 2, 3),)
    rules.append((q, atom_fix, bonds_fix, False))

    # Ozone
    #
    # [O*] -- O -- [O*]  >>  O == [O+] -- [O-]
    #
    q = smarts('[O;D1][O;D2][O;D1] |^1:0,2|')
    atom_fix = {1: (0, False), 2: (1, None), 3: (-1, False)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    # only after [A-] - [C+] rules!
    #  [C+] - N(R2)  >> C = [N+](R2)
    #
    q = smarts('[C;D1,D2,D3;z1;+]-[N;D3;z1;x0]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    # low priority for hetero-N
    #  [C+] - N(X2)  >> C = [N+](X2)
    #
    q = smarts('[C;D1,D2,D3;z1;+]-[N;D3;z1]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  [C+] = N(R)  >>  C # [N+](R)
    #
    q = smarts('[C;D1,D2;z2;x1;+]=[N;D1,D2;z2]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 3),)
    rules.append((q, atom_fix, bonds_fix, False))

    # fix rdkit
    # A-[Cl,Br,I+] - [O-] >> X = O
    #
    q = smarts('[Cl,Br,I;D2;z1;+][O;D1;-]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    # fix rdkit
    # A-[Hal+2]([O-])2
    #
    q = smarts('[Cl,Br,I;D3;z1;+2]([O;D1;-])[O;D1;-]')
    atom_fix = {1: (-2, None), 2: (1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    # fix rdkit
    # A-[Hal+3]([O-])3
    #
    q = smarts('[Cl,Br,I;D4;z1;+3]([O;D1;-])([O;D1;-])[O;D1;-]')
    atom_fix = {1: (-3, None), 2: (1, None), 3: (1, None), 4: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2), (1, 4, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    # fix reaxys [Cl-]=O > Cl-[O-]
    q = smarts('[Cl,Br,I;D1;z2;-]=[O;D1]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    return [(q, atom_fix, bonds_fix,
             [n for n, a in q.atoms() if a.atomic_symbol == 'A' and n not in atom_fix], is_tautomer)
            for q, atom_fix, bonds_fix, is_tautomer in rules]


rules = Proxy(_rules)


__all__ = ['rules']
