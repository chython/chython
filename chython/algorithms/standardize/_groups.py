# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    from ...containers import QueryContainer

    raw_rules = []

    #
    # A   H   A     A     H   A
    #  \ / \ /       \  .. \ /
    #   B   B    >>   B     B
    #  / \ / \       / \  .. \
    # A   H   A     A   H     A
    #
    atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'B', 'neighbors': 4, 'hybridization': 1},
             {'atom': 'H', 'neighbors': 2}, {'atom': 'H', 'neighbors': 2})
    bonds = ((1, 3, 1), (1, 4, 1), (2, 3, 1), (2, 4, 1))
    atom_fix = {}
    bonds_fix = ((1, 3, 8), (2, 4, 8))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #      A            A
    #     //           //
    # B - N >> [B-] - [N+]
    #     \            \
    #      A            A
    #
    atoms = ({'atom': 'B'}, {'atom': 'N', 'neighbors': 3, 'hybridization': 2})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #     A   A             A      A
    #     |   |             |      |
    # A - B = N - A >> A - [B-] - [N+] - A
    #     |   |             |      |
    #     A   A             A      A
    #
    atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 2}, {'atom': 'N', 'neighbors': 4, 'hybridization': 2})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #   A      A       A   A
    #   |      |       |   |
    #  [B-] = [N+] >>  B - N
    #   |      |       |   |
    #   A      A       A   A
    #
    atoms = ({'atom': 'B', 'charge': -1, 'neighbors': 3, 'hybridization': 2},
             {'atom': 'N', 'charge': 1, 'neighbors': 3, 'hybridization': 2})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #        [A-]                 A
    #         |                   |
    # [A-] - [B+3] - [A-] >> A - [B-] - A
    #         |                   |
    #        [A-]                 A
    #
    atoms = ({'atom': 'B', 'charge': 3, 'neighbors': 4}, {'atom': 'A', 'charge': -1}, {'atom': 'A', 'charge': -1},
             {'atom': 'A', 'charge': -1}, {'atom': 'A', 'charge': -1})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
    atom_fix = {1: (-4, None), 2: (1, None), 3: (1, None), 4: (1, None), 5: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #        A             A
    #        |             |
    # [A-] - B - A >> A - [B-] - A
    #        |             |
    #        A             A
    #
    atoms = ({'atom': 'B', 'neighbors': 4}, {'atom': 'A', 'charge': -1}, {'atom': 'A'}, {'atom': 'A'}, {'atom': 'A'})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #      A             A
    #      |             |
    # A -  B - A >> A - [B-] - A
    #      |             |
    #      A             A
    #
    atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1},)
    bonds = ()
    atom_fix = {1: (-1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #     A         A
    #     |         |
    #     N    >>  [N+]
    #   / | \     / | \
    #  A  A  A   A  A  A
    #
    atoms = ({'atom': 'N', 'neighbors': 4, 'hybridization': 1},)
    bonds = ()
    atom_fix = {1: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #  aci-nitro
    #      O          O
    #     //         //
    # C = N  >> C - [N+]
    #      \         \
    #       OH       [O-]
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1}, {'atom': 'C'})
    bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 2))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 4, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #            O               [O-]
    #           //               /
    # [C,N,O] = N  >> [C,N,O] = [N+]
    #           \                \
    #            A                A
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': ListElement(['C', 'N', 'O'])},
             {'atom': 'A'})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # fix aci-nitro
    #       [O-]          [O-]
    #       /             /
    # C = [N+ ] >>  C - [N+]
    #      \             \\
    #       OH            O
    #
    atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
             {'atom': 'O', 'neighbors': 1}, {'atom': 'C'})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 2))
    atom_fix = {}
    bonds_fix = ((1, 3, 2), (1, 4, 1))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # fix CN(=O)=N(=O)C
    #
    # [N+] = N = O >> [N+] = [N+] - O-
    #        |                |
    #        A                A
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'N', 'charge': 1, 'hybridization': 2, 'neighbors': 3}, {'atom': 'A'})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Fix CN(=O)=N(=N)C
    # [N+] = N = N - ? >> [N+] = [N+] - [N-] - ?
    #        |                    |
    #        A                    A
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
             {'atom': 'N', 'charge': 1, 'hybridization': 2, 'neighbors': 3}, {'atom': 'A'})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # For N-case not unique!
    #          N             [N-]
    #         //             /
    # [C,N] = N  >> [C,N] = [N+]
    #         \              \
    #          A              A
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
             {'atom': ListElement(['C', 'N'])}, {'atom': 'A'})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [N-] - [N+] = O >> N = [N+] - [O-]
    #         |               |
    #         A               A
    #
    atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3, 'hybridization': 2}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'N', 'charge': -1, 'hybridization': 1, 'neighbors': (1, 2)})
    bonds = ((1, 2, 2), (1, 3, 1))
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #         O              [O-]
    #        //              /
    # [A-] - N   >> [A-] - [N+]
    #        \\             \\
    #         O              O
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'A', 'charge': -1})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # O : N : O      O = [N+] - [O-]
    #     |      >>       |
    #     A               A
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1}, {'atom': 'A'})
    bonds = ((1, 2, 4), (1, 3, 4), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Nitrite
    #
    #   O        [O-]
    #  //        /
    # [N-]  >>  N
    #  \\       \\
    #   O        O
    #
    atoms = ({'atom': 'N', 'neighbors': 2, 'charge': -1}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 2), (1, 3, 2))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [O,C,N] = N # N >> [O,C,N] = [N+] = [N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': ListElement(['O', 'C', 'N'])})
    bonds = ((1, 2, 3), (1, 3, 2))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [C,N,O-] - [N+] # N  >> [C,N,O] = [N+] = [N-]
    #
    atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2}, {'atom': 'N', 'neighbors': 1},
             {'atom': ListElement(['N', 'C', 'O']), 'charge': -1, 'hybridization': 1})
    bonds = ((1, 2, 3), (1, 3, 1))
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  A - [N+] # N = [N-]  >> A - N = [N+] = [N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 2},
             {'atom': 'N', 'charge': -1, 'neighbors': 1}, {'atom': 'A'})
    bonds = ((1, 2, 3), (1, 3, 2), (2, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - N = N = N >> A - N = [N+] = [N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 2, 'hybridization': 2},
             {'atom': 'N', 'neighbors': 1})
    bonds = ((1, 2, 2), (1, 3, 2))
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - NH - N # N >> A - N = [N+] = [N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 2, 'hybridization': 1},
             {'atom': 'N', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 3))
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [N-] # N = N - A >> [N-] == [N+] == N - A
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1, 'charge': -1},
             {'atom': 'N', 'hybridization': 2, 'neighbors': 2})
    bonds = ((1, 2, 3), (1, 3, 2))
    atom_fix = {1: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [N-] == N # N >> [N-] == [N+] == [N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1, 'charge': -1}, {'atom': 'N', 'neighbors': 1})
    bonds = ((1, 2, 2), (1, 3, 3))
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 3, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - C # N = NH >> A - [CH] = [N+] = [N-]
    #
    atoms = ({'atom': 'C', 'neighbors': (1, 2)}, {'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1})
    bonds = ((1, 2, 3), (2, 3, 2))
    atom_fix = {2: (1, None), 3: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # note: order dependent
    # A - C # N = [O,N] >> A - C # [N+] - [O,N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': ListElement(['O', 'N']), 'hybridization': 2},
             {'atom': 'C', 'neighbors': (1, 2)})
    bonds = ((1, 2, 2), (1, 3, 3))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [NH2,OH,SH] - N # C - A >> [NH,O,S-] - [N+] # C
    #
    atoms = ({'atom': ListElement(['O', 'N', 'S']), 'neighbors': 1}, {'atom': 'N', 'neighbors': 2},
             {'atom': 'C', 'neighbors': (1, 2)})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - [NH] - N # C >> A - [N-] - [N+] # C
    #
    atoms = ({'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'neighbors': 2},
             {'atom': 'C', 'neighbors': (1, 2)})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # [NH2,OH,SH] - [N+] # [C-] >> [NH,O,S-] - [N+] # [CH]
    #
    atoms = ({'atom': ListElement(['O', 'N', 'S']), 'neighbors': 1}, {'atom': 'N', 'charge': 1, 'neighbors': 2},
             {'atom': 'C', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {1: (-1, None), 3: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - [NH] - [N+] # [C-] >> A - [N-] - [N+] # [CH]
    #
    atoms = ({'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'charge': 1, 'neighbors': 2},
             {'atom': 'C', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {1: (-1, None), 3: (1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # A - N # C >> A - [N+] # [C-]
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'C', 'neighbors': 1}, {'atom': 'A'})
    bonds = ((1, 2, 3), (1, 3, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # fix old diazo rule
    #
    # A - [C-] = [N+] = [NH] >> A - [CH] = [N+] = [N-]
    #
    atoms = ({'atom': 'C', 'charge': -1, 'neighbors': (1, 2), 'hybridization': 2},
             {'atom': 'N', 'charge': 1, 'neighbors': 2}, {'atom': 'N', 'neighbors': 1})
    bonds = ((1, 2, 2), (2, 3, 2))
    atom_fix = {1: (1, None), 3: (-1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #    |           |
    #  - N -   >> - [N+] -
    #    \\          |
    #    [O,N]      [O,N-]
    #
    atoms = ({'atom': 'N', 'neighbors': 4, 'hybridization': 2}, {'atom': ListElement(['O', 'N']), 'hybridization': 2})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # N-oxide radical
    #
    #    |         |
    #  - N*  >>  - N
    #    \\        |
    #     O        O*
    #
    atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 2, 'is_radical': True}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (0, False), 2: (0, True)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # C           C
    #  \           \
    #   N # N >>   [N+] = [N-]
    #  /           /
    # C           C
    #
    atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': 1}, {'atom': 'C'}, {'atom': 'C'})
    bonds = ((1, 2, 3), (1, 3, 1), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  C - N = [N+]  >>  C - [N+] # N
    #
    atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 1}, {'atom': 'C'})
    bonds = ((1, 2, 2), (1, 3, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 3),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  [N+] - [C-] = O  >>  N = C = O
    #
    atoms = ({'atom': 'C', 'charge': -1, 'neighbors': 2},
             {'atom': 'N', 'charge': 1, 'neighbors': (1, 2), 'hybridization': 1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  N # C - OH  >>  HN = C = O
    #
    atoms = ({'atom': 'C', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 3), (1, 3, 1))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  N # C - [O-]  >>  [N-] = C = O
    #
    atoms = ({'atom': 'C', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 3), (1, 3, 1))
    atom_fix = {2: (-1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # - [N+] - [O-]  >>  - N = O
    #
    atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2, 'hybridization': 1},
             {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  [C+] - N(R2)  >> C = [N+](R2)
    #
    atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 1},
             {'atom': 'C', 'charge': 1, 'hybridization': 1, 'neighbors': (1, 2, 3), 'heteroatoms': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  [C+] = N(R)  >>  C # [N+](R)
    #
    atoms = ({'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
             {'atom': 'C', 'charge': 1, 'neighbors': (1, 2), 'hybridization': 2, 'heteroatoms': 1})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 3),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #      [O,S]H    [O,S]
    #      /         //
    # N = C  >> NH - C
    #
    atoms = ({'atom': 'C', 'hybridization': 2}, {'atom': ListElement(['O', 'S']), 'neighbors': 1},
             {'atom': 'N', 'hybridization': 2})
    bonds = ((1, 2, 1), (1, 3, 2))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 1))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # fix pyridin-2-one. note: only after amide rule
    #
    atoms = ({'atom': 'C', 'neighbors': 3}, {'atom': ListElement(['O', 'S']), 'neighbors': 1},
             {'atom': 'N', 'neighbors': 2}, {'atom': 'A', 'hybridization': 2}, {'atom': 'A', 'hybridization': 2},
             {'atom': 'A', 'hybridization': 2}, {'atom': 'A', 'hybridization': 2})
    bonds = ((1, 2, 2), (1, 3, 1), (1, 4, 1), (3, 7, 1), (4, 5, (1, 2)), (5, 6, (1, 2)), (6, 7, (1, 2)))
    atom_fix = {}
    bonds_fix = ((1, 2, 1), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

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
    atoms = ({'atom': 'O', 'neighbors': 1}, {'atom': ListElement(['O', 'N'])}, {'atom': 'C'}, {'atom': 'C'})
    bonds = ((1, 3, 1), (2, 3, 1), (3, 4, 2))
    atom_fix = {}
    bonds_fix = ((1, 3, 2), (3, 4, 1))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #       A                   A
    #       |                   |
    #  A - [P+] - [O-]  >>  A - P = O
    #       |                   |
    #       A                   A
    #
    atoms = ({'atom': 'P', 'charge': 1, 'neighbors': 4, 'hybridization': 1},
             {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #       A                   A
    #       |                   |
    #  A - [P-] - [C+]  >>  A - P = C
    #       |                   |
    #       A                   A
    #
    atoms = ({'atom': 'P', 'charge': -1, 'neighbors': 4, 'hybridization': 1},
             {'atom': 'C', 'charge': 1, 'neighbors': (1, 2, 3), 'hybridization': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #   F   F         F   F
    #    \ /           \ /
    # F - P - F >> F - [P-] - F
    #    / \           / \
    #   F   F         F   F
    #
    atoms = ({'atom': 'P', 'neighbors': 6}, {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1},
             {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1},
             {'atom': 'F', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1))
    atom_fix = {1: (-1, None)}
    bonds_fix = ()
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #           A                   A
    #          /                   /
    # [O-] - [S,Si,Se+]  >>  O = [S,Si,Se]
    #          \                   \
    #           A                   A
    #
    atoms = ({'atom': ListElement(['S', 'Se', 'Si']), 'charge': 1, 'neighbors': 3, 'hybridization': 1},
             {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #       A                   A
    #       |                   |
    #  A = [S+] - [O-]  >>  A = S = O
    #       |                   |
    #       A                   A
    #
    atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 4, 'hybridization': 2},
             {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #       [O-]          O
    #       /            //
    # A = [S+2]  >>  A = S
    #       \            \\
    #       [O-]          O
    #
    atoms = ({'atom': 'S', 'charge': 2, 'neighbors': 3, 'hybridization': 2},
             {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 1))
    atom_fix = {1: (-2, None), 2: (1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #          A                  A
    #          |                  |
    # [O-] - [S+2] - [O-] >>  O = S = O
    #          |                  |
    #          A                  A
    #
    atoms = ({'atom': 'S', 'charge': 2, 'neighbors': 4, 'hybridization': 1},
             {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'O', 'charge': -1, 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 1))
    atom_fix = {1: (-2, None), 2: (1, None), 3: (1, None)}
    bonds_fix = ((1, 2, 2), (1, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  A - [S-] - [C+]  >>  A - S = C
    #       |                   |
    #       A                   A
    #
    atoms = ({'atom': 'S', 'charge': -1, 'neighbors': 3, 'hybridization': 1},
             {'atom': 'C', 'charge': 1, 'hybridization': 1, 'neighbors': (1, 2, 3)})
    bonds = ((1, 2, 1),)
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  O              O
    #  \\             \\
    #  [S-] - A  >>    S - A
    #  //             /
    #  O            [O-]
    #
    atoms = ({'atom': 'S', 'charge': -1, 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'A'})
    bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #  [O-] - [S+] = C  >>  O = S = C
    #
    atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 2}, {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'C'})
    bonds = ((1, 2, 1), (1, 3, 2))
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    #      A{1,3}               A{1,3}
    #      |                    |
    #  N = S - [OH]  >>  [NH] - S = O
    #
    atoms = ({'atom': 'S', 'neighbors': (3, 5), 'hybridization': 2}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2})
    bonds = ((1, 2, 1), (1, 3, 2))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 1))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # C # C - [O,NH,S]H  >> C=C=[O,NH,S]
    #
    atoms = ({'atom': ListElement(['O', 'S', 'N']), 'neighbors': 1}, {'atom': 'C', 'neighbors': 2},
             {'atom': 'C', 'neighbors': (1, 2)})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    #
    # C # C - [NH]R  >> C=C=NR
    #
    atoms = ({'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'C', 'neighbors': 2},
             {'atom': 'C', 'neighbors': (1, 2)})
    bonds = ((1, 2, 1), (2, 3, 3))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (2, 3, 2))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Carbon Monoxide
    #
    # [CX1] = O  >> [С-] # [O+]
    #
    atoms = ({'atom': 'C', 'neighbors': 1, 'is_radical': True}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 2),)
    atom_fix = {1: (-1, False), 2: (1, None)}
    bonds_fix = ((1, 2, 3),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Ozone
    #
    # [O*] -- O -- [O*]  >>  O == [O+] -- [O-]
    #
    atoms = ({'atom': 'O', 'neighbors': 1, 'is_radical': True}, {'atom': 'O', 'neighbors': 2},
             {'atom': 'O', 'neighbors': 1, 'is_radical': True})
    bonds = ((1, 2, 1), (2, 3, 1))
    atom_fix = {1: (0, False), 2: (1, None), 3: (-1, False)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Allyl radical
    #
    # [C*] -- [C*] -- [C*]  >>  C == C -- [C*]
    #
    atoms = ({'atom': 'C', 'is_radical': True}, {'atom': 'C', 'is_radical': True}, {'atom': 'C', 'is_radical': True})
    bonds = ((1, 2, 1), (2, 3, 1))
    atom_fix = {1: (0, False), 2: (0, False)}
    bonds_fix = ((1, 2, 2),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    # Br-ion + I-ion
    #
    #        A           A
    #        |           |
    # [Br-].[I+] >> Br - I
    #        |           |
    #        A           A
    #
    atoms = ({'atom': 'Br', 'charge': -1, 'neighbors': 0},
             {'atom': 'I', 'charge': 1, 'neighbors': 2, 'hybridization': 1},)
    bonds = ()
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 1),)
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    compiled_rules = []
    for atoms, bonds, atom_fix, bonds_fix in raw_rules:
        q = QueryContainer()
        for a in atoms:
            q.add_atom(**a)
        for n, m, b in bonds:
            q.add_bond(n, m, b)
        compiled_rules.append((q, atom_fix, bonds_fix))
    return compiled_rules


def _rules_double():
    from ...containers import QueryContainer

    raw_rules = []

    #
    #     [OH]                O
    #      |                 //
    #  N = S = A  >>  [NH] - S = A
    #      |                 |
    #      A                 A
    #
    atoms = ({'atom': 'S', 'neighbors': 4}, {'atom': 'O', 'neighbors': 1},
             {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2}, {'atom': 'A'}, {'atom': 'A'})
    bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 2), (1, 5, 1))
    atom_fix = {}
    bonds_fix = ((1, 2, 2), (1, 3, 1))
    raw_rules.append((atoms, bonds, atom_fix, bonds_fix))

    compiled_rules = []
    for atoms, bonds, atom_fix, bonds_fix in raw_rules:
        q = QueryContainer()
        for a in atoms:
            q.add_atom(**a)
        for n, m, b in bonds:
            q.add_bond(n, m, b)
        compiled_rules.append((q, atom_fix, bonds_fix))
    return compiled_rules


single_rules = Proxy(_rules_single)
double_rules = Proxy(_rules_double)


__all__ = ['single_rules', 'double_rules']
