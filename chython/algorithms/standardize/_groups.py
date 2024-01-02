# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...periodictable import ListElement


def _rules_single():
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
    #      [O,S]H    [O,S]
    #      /         //
    # N = C  >> NH - C
    #
    q = smarts('[N;z2]=[C;D2,D3;z2]-[O,S;D1]')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-2-one. note: only after amide rule
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[N;D2;z1][A;z2]-,=[A;z2][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-2-imine
    #
    q = smarts('[N;D2;z2;!R]=[C;D3;r6]1[N;D2;z1][A;z2]-,=[A;z2][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin-4-one
    #
    q = smarts('[O,S;D1;z2;x0]=[C;D3;r6]1[A;z2]=[A;z2][N;D2;z1][A;z2]-,=[A;z2]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    # todo:
    # [C;a:10][N;H:2][N:3]=[C:4]1[C:5]=,:[C:6][C:7](=[O:1])[C:8]=,:[C:9]1
    # [C;a:10][N;H:2][N:3]=[C:4]1[C:5](=[O:1])[C:6]=,:[C:7]-,:[C:8]=,:[C:9]1

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
    # fix pyridin. note: don't move.
    #
    q = smarts('[O,S,N;D1;z2;x0]=[C;D3;r6]1[N;D2;z2]=[A;z2][A;z2]-,=[A;z2][C;D2,D3;z1]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 7, 2))
    rules.append((q, atom_fix, bonds_fix, True))

    #
    # fix pyridin amine.
    #
    q = smarts('[N;D2;z2;!R]=[C;D3;r6]1[N;D2;z2]=[A;z2][A;z2]-,=[A;z2][C;D2,D3;z1]1')
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (2, 7, 2))
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


def _rules_double():
    from ... import smarts

    rules = []

    #
    #     [OH]                O
    #      |                 //
    #  N = S = A  >>  [NH] - S = A
    #      |                 |
    #      A                 A
    #
    q = smarts('[S;D4;z3:1]([O;D1:2])(=[N;D1,D2;z2:3])(=[A])[A]')
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 1))
    rules.append((q, atom_fix, bonds_fix, True))

    return [(q, atom_fix, bonds_fix,
             [n for n, a in q.atoms() if a.atomic_symbol == 'A' and n not in atom_fix], is_tautomer)
            for q, atom_fix, bonds_fix, is_tautomer in rules]


single_rules = Proxy(_rules_single)
double_rules = Proxy(_rules_double)


__all__ = ['single_rules', 'double_rules']
