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
"""The deliberately-nasty corpus the standardization pass is measured against.

79 public structures, every one drawn wrong on purpose.  It is the only asset reaching the 9
metal-organic rules; `test_standardize_groups_port.py` covers the functional groups and names the
unreached rules in its `UNREACHED`.  Kept as a SMILES literal rather than read from
`test/standardize.sdf`, so the gate depends on nothing but a parser while the format epic moves.
"""

CORPUS = (
    '[H]1B(C)([H]B1(C)C)C',
    'B([N](=C)C)(C)(C)C',
    'N(C)(=B(C)(C)C)(C)C',
    '[N](C)(B(C)(C)C)(C)C',
    '[N+](C)(=[B-](C)C)C',
    '[B+3]([F-])([F-])([F-])[F-]',
    'B([F-])(F)(F)F',
    'O([O])[O] |^1:1,2|',
    'C=N(=O)O',
    'C=N(C)=O',
    'CN(=N)=O',
    'CN(=N)=N',
    'C=N(C)=N',
    'C[N+]([NH-])=O',
    'CN(=O)=N(C)=N',
    'O=N(C)=N(=O)C',
    'C=[N+]([O-])O',
    'CN(=O)=O',
    'N([O-])(=O)=O',
    'Cn(o)o',
    '[N-](=O)=O',
    'C=N#N',
    'N#N=O',
    'N#N=N',
    '[CH2-][N+]#N',
    'N#[N+][O-]',
    '[NH-][N+]#N',
    'C[N+]#N=[N-]',
    'CN=[N]=N',
    'CN[N]#N',
    'CN=N#[N-]',
    '[N-]=N#N',
    'C#N=N',
    'C#N=NC',
    'C#N=O',
    '[CH-]=[N+]=N',
    'CN(=N)(C)C',
    'CN(=O)(C)C',
    'C[N](=O)C |^1:1|',
    'CN(#N)C',
    'CN=[NH2+]',
    '[C-]([NH2+]C)=O',
    'C(#N)O',
    'C(#N)[O-]',
    'C[NH2+][O-]',
    '[CH+](C)N(C)C',
    '[C+](=N\\C)/C',
    'C(=N\\C)/O',
    'C[P+]([O-])(C)C',
    'C[P-]([CH2+])(C)C',
    'FP(F)(F)(F)(F)F',
    'C[S+]([O-])C',
    'C=[S+](C)([O-])C',
    'C=[S+2]([O-])[O-]',
    'C[S+2]([O-])(C)[O-]',
    'C[S-]([CH2+])C',
    'C[S-](=O)=O',
    'C=[S+][O-]',
    'CS(=N)O',
    'CN=S(=N)(O)O',
    '[CH]=O |^1:0|',
    'C(#C)O',
    'C(#C)NC',
    '[Br-].Br[I+]Br',
    'C#[N]O',
    'C#[N]NC',
    '[C-]#[N+]NC',
    '[C-]#[N+]O',
    'C#[N]OC',
    'C=1N(C)C(=[Cu]I)N(C)C=1',
    'C(#N)[Fe]',
    '[C-](#N)[Fe]',
    'C(=O)[Fe](C=O)C=O',
    'C1(=O)[Fe]C(=O)[Fe]1',
    '[C](=O)[Fe] |^1:0|',
    'C(#O)[Fe]',
    '[H]C1=2C=3([H])[Fe+2]4156789(C=1([H])C6(=C4([H])[C-]5([H])C=17[H])[H])C=3([H])[C-]9([H])C=28[H]',
    'C12=C3[Fe+2]1456789(C=1([Li])C5=C8[C-]4C=19)C(=C36)[C-]27',
    '[Fe]12345678(C9C1C6C4C39)C1C2C7C5C81',
)


__all__ = ['CORPUS']
